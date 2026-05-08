"""
Estimate complexity and run full key recovery attack when feasible.

usage: attack-aim2.py [-h] [--pairs-file PAIRS_FILE] [-n N] [-l L] [-L L] [--only-estimate] [--seed SEED] [--no-progress-bars]

options:
  -h, --help            show this help message and exit
  --pairs-file PAIRS_FILE
                        Path to file with (iv,ct) pairs
  -n N                  field dimension n (bits), default: deduced from file, or 128 if random
  -l L, -ell L          number of branches, default: 2 if n in {128,192}, 3 if n == 256
  -L L                  number of IVs, default: all if from file, or 1 if random
  --only-estimate       Do not run the full attack, only show complexity estimates
  --seed SEED           Randomness seed
  --no-progress-bars    Disable progress bar (for cleaner logs)


Option 1: generate random IVs to construct the system. ell is determined automatically, L=1 by default.

$ sage attack-aim2.py -n 128

$ sage attack-aim2.py -n 128 -L 50  # few seconds

$ sage attack-aim2.py -n 128 -L 15  # 1 hour

$ sage attack-aim2.py -n 128 -l 5  # override with 5 branches


Option 2: run on file with (iv, ct) pairs generated from C reference implementation (see keypairs_from_reference_implementation/simple_keypairs.c).

$ sage attack-aim2.py --pairs-file keypairs_from_reference_implementation/keypairs128.txt

$ sage attack-aim2.py --pairs-file keypairs_from_reference_implementation/keypairs256.txt

"""

import os, sys, psutil, math, time, pathlib, random, argparse, gc
from collections import namedtuple

from sage.all import binomial, RR, log, matrix, prod

try:
    from tqdm import tqdm
except ImportError:
    tqdm = lambda v: v

import ff_rust

from polynomial import PolynomialRing as PolynomialRing_Generic
from halfgcd import PolyGCD, PolyHalfGCD, PolyQuoRem, set_FAFFT

from aim2 import AIM2_Rust, parse_keypairs_file


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs-file", type=pathlib.Path, help="Path to file with (iv,ct) pairs")
    parser.add_argument("-n", type=int, help="field dimension n (bits), default: deduced from file, or 128 if random")
    parser.add_argument("-l", "-ell", type=int, help="number of branches, default: 2 if n in {128,192}, 3 if n == 256")
    parser.add_argument("-L", type=int, help="number of IVs, default: all if from file, or 1 if random")
    parser.add_argument("--related-key", action="store_true", help="Instead of changing the IVs, compute on related keys.")
    parser.add_argument("--only-estimate", action="store_true", help="Do not run the full attack, only show complexity estimates")
    # parser.add_argument("--pre-reduce", action="store_true", help="Pre-reduce matrix")
    parser.add_argument("--seed", type=int, help="Randomness seed")
    parser.add_argument("--no-progress-bars", action="store_true", help="Disable progress bar (for cleaner logs)")
    args = parser.parse_args()

    if args.seed is None:
        args.seed = int(os.urandom(8).hex(), 16)

    if args.no_progress_bars:
        global tqdm
        tqdm = lambda v: v

    print("[i] Seed:", args.seed)
    random.seed(args.seed)

    if args.pairs_file:
        if args.n is not None or args.l is not None:
            print("ERROR: n or l can not be specified for pairs file")
            print()
            parser.print_usage()
            return
        if args.related_key:
            print("ERROR: related key attack on pairs file not supported")
            print()
            parser.print_usage()
            return

        print("Reading secret key (optional) and iv/ct pairs from file", args.pairs_file, "...")
        SECRET_KEY, iv_ct_pairs = parse_keypairs_file(args.pairs_file)
        args.n = len(SECRET_KEY) * 8
        args.l = 2 if args.n < 224 else 3
        if args.L is None:
            args.L = len(iv_ct_pairs)
        assert args.L <= len(iv_ct_pairs), "not enough pairs"

        AIM = AIM2_Rust.BY_N[args.n]
        SECRET_KEY = AIM.gf_from_bytes(SECRET_KEY)
        print("[i] Secret key:", repr(SECRET_KEY))
        print("\nVerifying secret key versus %d (iv, ct) pairs..." % len(iv_ct_pairs))
        for i, (iv, ct) in enumerate(tqdm(iv_ct_pairs)):
            ct2 = AIM(iv).eval(SECRET_KEY)
            if ct != ct2:
                print("[W] Secret key verification failed on pair #%d, stopping checks..." % i)
                break
        else:
            print("[i] Secret key supplied for verification matches the (iv, ct) pairs")
        print()

    else:
        if args.n is None:
            args.n = 128
        if args.l is None:
            args.l = 2 if args.n < 224 else 3
        if args.L is None:
            args.L = 1

        if args.related_key:
            assert args.l == 2, "Related-key only supported for ell=2"
            assert args.L <= args.n

        AIM = AIM2_Rust.custom_instance(n=args.n, ell=args.l)
        SECRET_KEY = AIM.F(random.randrange(2**AIM.n))
        iv = random.randrange(2**AIM.n).to_bytes(AIM.nb)  # for related-key only
        A = AIM(iv=iv)
        iv_ct_pairs = []
        print("[i] Secret key:", repr(SECRET_KEY))
        print("\nGenerating (iv, ct) pairs...")
        for i in tqdm(range(args.L)):
            if args.related_key:
                γ1, γ2 = AIM.cs
                a = AIM.F(1) / (γ1 + γ2)
                b = AIM.cs[0] * a

                sk_i = ((a * SECRET_KEY + b)**(2**i) + b) / a
                ct_i = A.eval(sk_i)

                iv_ct_pairs.append((iv, (ct_i, i)))
            else:
                iv = random.randrange(2**AIM.n).to_bytes(AIM.nb)
                ct = AIM(iv).eval(SECRET_KEY)
                iv_ct_pairs.append((iv, ct))
        print()

    run_attack(
        args.n, args.l, args.L, AIM, iv_ct_pairs,
        related_key=args.related_key,
        only_estimate=args.only_estimate,
    )
    print("[i] Real secret key for comparison:", repr(SECRET_KEY))


def create_initial_equations(AIM, iv_ct_pairs, p, ts):
    eqs_base = [
        (p + AIM.cs[i]) * ts[i] - ts[i]**(2**AIM.es[i])
        for i in range(AIM.ell)
    ]
    print("\nGenerating base equations... (using linearized polynomials)")
    for iv, ct in tqdm(iv_ct_pairs):
        ct = AIM.gf_from_bytes(ct)
        A = AIM(iv)

        Lpoly = A.LIN_const
        for branch in range(A.ell):
            Lpoly += A.apply_linearized(ts[branch], A.LIN_polys[branch])

        fl = Lpoly*(p + ct) - Lpoly**(2**AIM.es[-1])
        eqs_base.append(fl)
    print()
    return eqs_base


def create_initial_equations_related_key(AIM, iv_ct_pairs, p, ts):
    assert AIM.ell == 2
    γ1, γ2 = AIM.cs
    a = AIM.F(1) / (γ1 + γ2)
    b = AIM.cs[0] * a

    assert a * AIM.cs[0] + b == 0
    assert a * AIM.cs[1] + b == 1

    eqs_base = [
        (p + 0) * ts[0] - ts[0]**(2**AIM.es[0]),
        (p + 1) * ts[1] - ts[1]**(2**AIM.es[1]),
    ]

    print("\nGenerating base RK equations... (using linearized polynomials)")
    for iv, (ct_i, i) in tqdm(iv_ct_pairs):
        ct_i = AIM.gf_from_bytes(ct_i)
        new_ct = (a * ct_i + b)**(2**(AIM.n - i))

        A2s = AIM(iv=iv).with_affine_shift(a, b).with_squaring(i)

        Lpoly = A2s.LIN_const
        Lpoly += A2s.apply_linearized(ts[0], A2s.LIN_polys[0])
        Lpoly += A2s.apply_linearized(ts[1], A2s.LIN_polys[1])

        fl = Lpoly*(p + new_ct) - Lpoly**(2**AIM.es[-1])
        eqs_base.append(fl)

        # A2s = A2.with_squaring(i)
        # ct = A2s.eval(sk, as_bytes=False)
        # test = (a * ct_i + b)**(2**(A.n-i))
    print()
    return eqs_base



def square_equations(eqs, M):
    eqs2 = []
    for eq in eqs:
        for _ in range(M+1):
            eqs2.append(eq)
            eq = eq**2
    return eqs2


def collect_t_monomials(eqs, p):
    t_monomials = set()
    for eq in eqs:
        for _, mono in eq:
            mono_t = mono * p**(-mono.degrees()[0])
            t_monomials.add(mono_t)
    t_monomials = sorted(t_monomials, reverse=True)
    return {mono: i for i, mono in enumerate(t_monomials)}


RowInfo = namedtuple("RowInfo", ("degree", "nnz", "eq_id", "row_const", "row_full"))


def equations_to_polynomial_matrix(eqs, Rx, t_monomial_index, ts, p0, only_estimate):
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


def run_attack(n, ell, L, AIM, iv_ct_pairs, related_key=False, only_estimate=False):
    M = math.ceil( (ell * n + 1) / (L + ell)) - 1
    W = n*ell + 1

    print("\n[i] AIM/attack parameters:")
    print(f"- n = {n} bits")
    print(f"- ell = {ell} branches")
    print(f"- L = {L} IVs")
    print(f"- M = {M} squarings (1 ... 2^(M-1) 2^M)")
    print(f"- W = {W} t-monomials")
    print(f"- Only estimate ? {only_estimate}")
    # print(f"- Pre-reduce ? {pre_reduce}")
    print()

    time_comp_D0 = 21 * 2**(ell * n / (ell + L)) * (ell + L) * n**3
    print("[I] Time complexity rough estimate = 2^%.2f field ops" % math.log(time_comp_D0, 2))
    assert n >= 8
    assert ell >= 1

    #F = AIM.F
    F = getattr(ff_rust, "GF%dField" % n)()
    Rx = getattr(ff_rust, "PolynomialRingGF%d" % n)()
    FAcls = getattr(ff_rust, "FAFFT_GF%d" % n)
    PivotRows = getattr(ff_rust, "pivot_rows_gf%d" % n)
    Inverse = getattr(ff_rust, "inverse_dense_gf%d" % n)
    SingleApproximant = getattr(ff_rust, "single_approximant_gf%d" % n)

    names = "p"
    for i in range(AIM.ell):
        names += ",t%d" % (i+1)
    p, *ts = PolynomialRing_Generic(F, names=names).gens()

    if related_key:
        eqs = create_initial_equations_related_key(AIM, iv_ct_pairs, p, ts)
    else:
        eqs = create_initial_equations(AIM, iv_ct_pairs, p, ts)
    print("[i] Initial equations:", len(eqs))
    eqs = square_equations(eqs, M)
    print("[i] Squared equations:", len(eqs))
    H = len(eqs)
    assert H == len(eqs) >= W, "not enough equations.."

    t_monomial_index = collect_t_monomials(eqs, p)
    print("[i] W = t-monomials in equations:", len(t_monomial_index))
    assert W == len(t_monomial_index), "mispredicted #monomials.."

    eqs = iter(eqs)  # free after consumption
    row_infos = equations_to_polynomial_matrix(eqs, Rx, t_monomial_index, ts, p0=F(0), only_estimate=only_estimate)
    row_infos.sort(key=lambda ri: (ri.degree, ri.nnz, ri.eq_id))
    assert len(row_infos) == H

    print("[i] Full equation matrix F: H=%d x W=%d" % (H, W))

    print("\nComputing pivots (p=0)")
    pivots = PivotRows([ri.row_const for ri in row_infos])

    print("[i] Rank = #pivots = %d/%d" % (len(pivots), W), "  fingerprint %x" % abs(hash(str(pivots))))
    assert len(pivots) == W, "non-full rank!"

    ksi = sum(row_infos[i].degree for i in pivots)
    ksi_estimated = (ell * n + 1) * 2**(M+1) / (M+1)

    print("[I] Row-degree sum ξ: experimental = 2^%.2f, estimated = 2^%.2f" % (
        math.log(ksi, 2),
        math.log(ksi_estimated, 2),
    ))
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
        f"[R1] n={n:3d}  ell={ell:2d}  L={L:3d}  M={M:3d}  W={W:4d}"
        f"  ξ=2^{math.log(ksi, 2):6.2f}"
        f"  Time=2^{math.log(time_comp, 2):6.2f} field ops"
        f"  Mem=2^{math.log(mem_comp, 2):6.2f} field elems"
    )
    print(
        f"[R2] n={n:3d}  ell={ell:2d}  L={L:3d}  M={M:3d}  W={W:4d}"
        f"  ξ=2^{math.log(ksi, 2):6.2f}"
        f"  Time=2^{math.log(time_comp/n, 2):6.2f} enc"
        f"  Mem=2^{math.log(mem_bytes, 2):6.2f} bytes"
    )
    print()

    avail_bytes = psutil.virtual_memory().available
    if only_estimate or mem_bytes >= 0.9*avail_bytes:
        print("[W] not running full attack", "%.1f GB needed, %.1f GB available" % (mem_bytes / 2**30.0, avail_bytes / 2**30.0))
        return


    print("[i] Running full attack", "%.1f GB needed, %.1f GB available" % (mem_bytes / 2**30.0, avail_bytes / 2**30.0))
    print()

    matF = [row_infos[i].row_full for i in pivots]

    # polynomials
    Fconst = [row[-1] for row in matF]
    Fplus = [row[:-1] for row in matF]

    # mod p (constants)
    Fplus0 = [[v[0] for v in row] for row in Fplus]


    # Step 1: Kernel vector
    G = Fplus
    G0 = Fplus0

    o = ksi * 2 + 2
    print("[o] Approximant order o = %d = 2^%.2f" % (o, math.log(o, 2)))

    FA_deg = int(math.ceil(math.log(o * 2 + 1, 2)))
    assert 2**FA_deg >= o * 2
    print("[i] FAFFT max degree 2^%d" % FA_deg)
    FA = FAcls(2**FA_deg)
    set_FAFFT(FA)

    x = Rx.x()
    mod = x**o

    print("[T] Main attack part, start time measurement", memory_usage())
    t0 = time.time()

    print("\n[.] Computing pivots of G0 for inversion...")
    pivots = PivotRows(G0)
    assert len(pivots) == W-1

    G0base = [G0[i] for i in pivots]

    print("[.] Inversion of base matrix...")
    G0base_inv = Inverse(G0base)

    i = row = None
    for i, row in enumerate(G0):
        if i not in pivots:
            break
    inds = list(pivots) + [i]

    print("[.] Single approximant...")
    #z0 = -row * G0_inv  # z * G0 = row
    z0 = [
        sum([row[i] * G0base_inv[i][j] for i in range(W-1)], F(0))
        for j in range(W-1)
    ]
    Gperm = [G[i] for i in inds]
    ker = SingleApproximant(Gperm, G0base_inv, o, z0)
    assert ker[-1] == Rx.constant(F(1))
    gc.collect()

    print("[T] TIME after single approximant:", time.time() - t0, memory_usage())
    print()

    print("[i] Kernel len =", len(ker), "degrees =", [v.degree() for v in ker])

    f = ker[0]  # heuristic
    #f = sum([AIM.F.random_element() * pol for pol in ker[:-1]], Rx.zero())  # probabilistic

    denom = PolyHalfGCD(mod, f, Rx=Rx, destroy=True)[-1]
    del mod, f
    gc.collect()
    print("[T] TIME after first HalfGCD:", time.time() - t0, memory_usage())
    degs_kerF = []

    Fconst_perm = [Fconst[i] for i in inds]
    pol = Rx.zero()
    for i in range(len(ker)):
        dv = FA.multiply_auto(denom, ker[i])
        ker[i] = None # free memory
        dv.resize(o)
        dv.in_place_trim()
        degs_kerF.append(dv.degree())
        pol += Fconst_perm[i] * dv
    del dv
    gc.collect()

    print("[i] Degrees of kernel * F", degs_kerF)

    print("[i] Got univariate poly in p, degree", pol.degree())

    print("[T] TIME to get univariate polynomial:", time.time() - t0, memory_usage())
    print()

    # manual roots... M(n)(log n)...
    # x^(2^n) - x  modulo  pol(x)
    cur = Rx.x()
    for i in range(n):
        cur = cur.square()
        if cur.degree() >= pol.degree():
            _, cur = PolyQuoRem(cur, pol, Rx=Rx)
    cur -= Rx.x()

    print("[T] TIME pre gcd (x^(2^n)-x mod f):", time.time() - t0, memory_usage())

    gpol = PolyGCD(pol, cur, Rx=Rx, destroy=True)
    del pol, cur
    print("gcd polynomial", gpol)
    print("[T] TIME post gcd:", time.time() - t0, memory_usage())

    gpol_sage = AIM.Fsage['x']([AIM.Fsage.from_integer(v.to_integer()) for v in gpol])
    print("found roots:")
    roots = gpol_sage.roots()
    for v, _ in roots:
        vs = ("%x" % v.to_integer()).zfill(n // 4)
        print("-", vs, end=" | ")

        v = F(v.to_integer())
        if not related_key:
            n_match = sum(AIM(iv=iv).eval(v) == ct for iv, ct in iv_ct_pairs)
            print(n_match, "/", len(iv_ct_pairs), "match")
        else:
            γ1, γ2 = AIM.cs
            a = AIM.F(1) / (γ1 + γ2)
            b = AIM.cs[0] * a
            vreal = (v - b) / a
            print("unshifted pt", repr(vreal))
    print("[T] TIME full attack (main part only):", time.time() - t0, memory_usage())


def memory_usage():
    used = psutil.Process().memory_info().rss
    used = max(1, used)
    return " | memory usage: %.2f GB" % (used/2**30.0)


if __name__ == '__main__':
    main()
