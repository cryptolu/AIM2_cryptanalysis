"""
Estimate complexity of attacks on RAIN.

usage: attack-rain.py [-h] [-n N] [-r R] [-L L] [-M M] [--seed SEED] [--no-progress-bars]

options:
  -h, --help          show this help message and exit
  -n N                field dimension n (bits), default: deduced from file, or 128 if random
  -r R                number of rounds
  -L L                number of IVs (unused for now)
  -M M                parameter M
  --seed SEED         Randomness seed
  --no-progress-bars  Disable progress bar (for cleaner logs)


2 rounds: (M automatic)

$ sage attack-rain.py -n 128 -r 2

$ sage attack-rain.py -n 192 -r 2

$ sage attack-rain.py -n 256 -r 2


3 rounds: (RAM consuming!!!)

$ sage attack-rain.py -n 32 -r 3 -M 23  # Toy version (32 bits)

$ sage attack-rain.py -n 128 -r 3 -M 109   # 1.5 hours 42GB RAM !

"""

import os, sys, psutil, math, time, pathlib, random, argparse, gc
from collections import namedtuple

from sage.all import binomial, RR, log, matrix, prod, ZZ, GF

try:
    from tqdm import tqdm
except ImportError:
    tqdm = lambda v: v

import ff_rust

from polynomial import PolynomialRing as PolynomialRing_Generic
from halfgcd import PolyGCD, PolyHalfGCD, PolyQuoRem, set_FAFFT

from aim2 import AIM2_Rust


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, default=128, help="field dimension n (bits), default: deduced from file, or 128 if random")
    parser.add_argument("-r", type=int, default=2, help="number of rounds")
    parser.add_argument("-L", type=int, default=1, help="number of IVs")
    parser.add_argument("-M", type=int, help="parameter M")
    # parser.add_argument("--only-estimate", action="store_true", help="Do not run the full attack, only show complexity estimates")
    parser.add_argument("--seed", type=int, help="Randomness seed")
    parser.add_argument("--no-progress-bars", action="store_true", help="Disable progress bar (for cleaner logs)")
    args = parser.parse_args()

    if args.seed is None:
        args.seed = int(os.urandom(8).hex(), 16)

    if args.no_progress_bars:
        global tqdm
        tqdm = lambda v: v

    if args.M is None:
        if args.r == 2:
            args.M = args.n//2 - 1

    # args.only_estimate = True

    print("[i] Seed:", args.seed)
    random.seed(args.seed)

    # use AIM to generate linear layers
    AIM = AIM2_Rust.custom_instance(n=args.n, ell=args.r)
    A = AIM(iv=b"")
    #Ls = AIM.LIN_polys
    Ls = []
    for mat in A.LIN_mats:
        mat = [ZZ(v.to_integer()).digits(2, padto=AIM.n) for v in mat]
        mat = matrix(GF(2), mat)
        imat = ~mat
        mat = [ AIM.Fsage(vec) for vec in mat]
        imat = [ AIM.Fsage(vec) for vec in imat]
        Lmat = AIM.matrix_to_linearized(mat)
        Limat = AIM.matrix_to_linearized(imat)
        Ls.append((Lmat, Limat))
    cs = [AIM.F(random.randrange(2**AIM.n)) for _ in range(args.r)]
    k = SECRET_KEY = AIM.F(random.randrange(2**AIM.n))
    assert AIM.apply_linearized(AIM.apply_linearized(k, Ls[0][0]), Ls[0][1]) == k

    iv_ct_pairs = []
    print("[i] Secret key:", repr(SECRET_KEY))
    print("\nGenerating (iv=c2, ct) pairs...")
    v = x = AIM.F(random.randrange(2**AIM.n))
    log = []
    for i in range(args.r-1):
        v += k + cs[i]
        v = v**(2**AIM.n-2)
        ti = v
        v = AIM.apply_linearized(v, Ls[i][0])
        t = v
        log.append((ti, t))
    v += k + cs[args.r-1]
    v = v**(2**AIM.n-2)
    v += k
    y = v
    iv_ct_pairs.append((x, y, log))

    # for _ in tqdm(range(args.L)):
    #     iv = random.randrange(2**AIM.n).to_bytes(AIM.nb)
    #     ct = AIM(iv).eval(SECRET_KEY)
    #     iv_ct_pairs.append((iv, ct))
    # print()

    run_attack(
        args.n, args.r, args.L, args.M, Ls, cs, iv_ct_pairs, AIM,
        # only_estimate=args.only_estimate,
    )
    print("[i] Real secret key for comparison:", repr(SECRET_KEY))


def square_equations(eqs, M):
    eqs2 = []
    for eq in eqs:
        for _ in range(M+1):
            eqs2.append(eq)
            eq = eq**2
    return eqs2


def collect_t_monomials(eqs, p):
    t_monomials = set()
    for eq in tqdm(eqs):
        for _, mono in eq:
            mono_t = mono * p**(-mono.degrees()[0])
            t_monomials.add(mono_t)
    t_monomials = sorted(t_monomials, reverse=True)
    return {mono: i for i, mono in enumerate(t_monomials)}


RowInfo = namedtuple("RowInfo", ("degree", "nnz", "eq_id", "row_const", "row_full"))


def equations_to_polynomial_matrix(eqs, Rx, t_monomial_index, ts, p0, only_estimate=True):
    x = Rx.x()
    W = len(t_monomial_index)
    F = type(p0)

    print("\nArranging equations...")
    rows = []
    for ieq, eq in enumerate(tqdm(eqs)):
        if not only_estimate:
            row_full = [Rx.zero() for _ in range(W)]
        else:
            row_full = None
        row_const = [F(0) for _ in range(W)]
        row_deg = -1
        nnz = 0
        for coef, mono in eq:
            deg_p, *deg_ts = mono.degrees()
            row_deg = max(row_deg, deg_p)
            mono_t = prod(ti**degi for ti, degi in zip(ts, deg_ts)) # mono * p**(-mono.degrees()[0])
            imono = t_monomial_index[mono_t]
            if coef:
                if not only_estimate:
                    row_full[imono] += coef * x**deg_p
                nnz += 1
                row_const[imono] += coef * p0**deg_p

        rows.append(RowInfo(degree=row_deg, nnz=nnz, eq_id=ieq, row_const=row_const, row_full=row_full))
    print()
    return rows


def run_attack(n, nr, L, M, Ls, cs, iv_ct_pairs, AIM, only_estimate=True):
    # note: AIM is only used to reuse linearized polynomial functions!
    print("\n[i] RAIN/attack parameters:")
    print(f"- n = {n} bits")
    print(f"- r = {nr} rounds")
    print(f"- L = {L} IVs")
    print(f"- M = {M} squarings (1 ... 2^(M-1) 2^M)")
    print(f"- Only estimate ? {only_estimate}")
    # print(f"- Pre-reduce ? {pre_reduce}")
    print()

    assert n >= 8
    assert nr >= 2

    #F = AIM.F
    F = getattr(ff_rust, "GF%dField" % n)()
    Rx = getattr(ff_rust, "PolynomialRingGF%d" % n)()
    PivotRows = getattr(ff_rust, "pivot_rows_gf%d" % n)

    names = "k"
    for i in range(nr-1):
        names += ",t%d" % (i+1)
    k, *ts = PolynomialRing_Generic(F, names=names).gens()

    x, y, _log = iv_ct_pairs[0]

    if nr == 2:
        t1, = ts
        c1, c2 = cs

        t1i = AIM.apply_linearized(t1, Ls[0][1])

        f1 = (x + k + c1) * t1i + AIM.F(1)
        f2 = (t1 + k + c2) * (k + y) + AIM.F(1)

        eqs = []
        eqs.append(f1)
        eqs.append(f1*t1i)

        eqs.append(f2)
        eqs.append(f2*t1)

        print("[i] Initial equations:", len(eqs))
        eqs = square_equations(eqs, M)
        print("[i] Squared equations:", len(eqs))
    else:
        t1, t2i = ts
        c1, c2, c3 = cs

        t1i = AIM.apply_linearized(t1, Ls[0][1])  # L1^inv
        t2 = AIM.apply_linearized(t2i, Ls[1][0])  # L2

        f1 = (x + k + c1) * t1i + AIM.F(1)
        f2 = (t1 + k + c2) * t2i + AIM.F(1)
        f3 = (t2 + k + c3) * (k + y) + AIM.F(1)

        t2i_sqs = [t2i]
        for _ in range(M):
            t2i_sqs.append(t2i_sqs[-1]**2)
        t1_sqs = [t1]
        for _ in range(M):
            t1_sqs.append(t1_sqs[-1]**2)

        base = [
            f1,
            f1*t1i,

            f2,
            f2*t1,
            f2*t2i,

            # f3 - not helpful?
            f3*t2,
        ]
        eqs = []
        for m in tqdm(range(M+1)):
            eqs.extend(base)
            for j in range(M+1):
                eqs.append(base[1] * t2i_sqs[j])
                eqs.append(base[-1] * t1_sqs[j])
                if m == 0:
                    eqs.append(f1 * t2i_sqs[j])
                    eqs.append(f3 * t1_sqs[j])
            if m != M:
                base = [f**2 for f in base]

    print("[i] Final equations:", len(eqs))
    H = len(eqs)

    t_monomial_index = collect_t_monomials(eqs, k)
    # for v in sorted(t_monomial_index):
    #     print(v._degrees)
    W = len(t_monomial_index)
    print("[i] W = t-monomials in equations:", len(t_monomial_index))

    assert H == len(eqs) >= W, "not enough equations.."

    assert W == len(t_monomial_index), "mispredicted #monomials.."

    #eqs = iter(eqs)  # free after consumption
    # p0 = F.random_element()
    p0 = F(0)
    row_infos = equations_to_polynomial_matrix(eqs, Rx, t_monomial_index, ts, p0=p0, only_estimate=only_estimate)
    row_infos.sort(key=lambda ri: (ri.degree, ri.nnz, ri.eq_id))
    assert len(row_infos) == H

    print("[i] Full equation matrix F: H=%d x W=%d nnz=%d" % (H, W, sum(ri.nnz for ri in row_infos)))

    print("\nComputing pivots (p=0)")
    pivots = PivotRows([ri.row_const for ri in row_infos])

    print("[i] Rank = #pivots = %d/%d" % (len(pivots), W), "  fingerprint %x" % abs(hash(str(pivots))))
    assert len(pivots) == W, "non-full rank!"

    ksi = sum(row_infos[i].degree for i in pivots)
    print("[I] Row-degree sum ξ: experimental = 2^%.2f" % math.log(ksi, 2))
    total_nnz = sum(row_infos[i].nnz for i in pivots)
    print("[I] Total nonzero coefficients", total_nnz, "<= 2W^2 =", 2*W**2, "?", total_nnz <= 2*W**2)

    Mksi = ksi * math.log(ksi, 2)**2 / 2
    time_comp = ksi * (4*W**2 + total_nnz) + Mksi * (2*n + 2*W + 22 * log(ksi, 2))

    print("[I] Time complexity 2^%.2f field ops = 2^%.2f encs  (1/n convention)" % (
        math.log(time_comp, 2), math.log(time_comp/n, 2),
    ))

    mem_comp = ksi * 2 * W + 10 * W**2
    mem_bytes = mem_comp * n / 8
    print("[I] Memory complexity : 2^%.2f field elems = 2^%.2f bytes = %.2f GB" % (
        math.log(mem_comp, 2), math.log(mem_bytes, 2), mem_bytes / 2**30.0
    ))
    print(
        f"[R1] n={n:3d}  r={nr:2d}  L={L:3d}  M={M:3d}  W={W:4d}"
        f"  ξ=2^{math.log(ksi, 2):6.2f}"
        f"  Time=2^{math.log(time_comp, 2):6.2f} field ops"
        f"  Mem=2^{math.log(mem_comp, 2):6.2f} field elems"
    )
    print(
        f"[R2] n={n:3d}  r={nr:2d}  L={L:3d}  M={M:3d}  W={W:4d}"
        f"  ξ=2^{math.log(ksi, 2):6.2f}"
        f"  Time=2^{math.log(time_comp/n, 2):6.2f} enc"
        f"  Mem=2^{math.log(mem_bytes, 2):6.2f} bytes"
    )
    print()


def memory_usage():
    used = psutil.Process().memory_info().rss
    used = max(1, used)
    return " | memory usage: %.2f GB" % (used/2**30.0)


if __name__ == '__main__':
    main()
