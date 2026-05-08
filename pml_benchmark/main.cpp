#include <iostream>
#include <algorithm>
#include <fstream>
#include <unordered_set>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/ZZ_pX.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/vec_ZZ_p.h>
#include <chrono>


/*
W 2 DEG 65536 TIME 4.5s
W 8 DEG 8192 TIME 2.7s
W 16 DEG 4096 TIME 3.7s
W 32 DEG 2048 TIME 5s
W 64 DEG 512 TIME 6s
W 128 DEG 64 TIME 4.7s
W 256 DEG 16 TIME 7s
W 512 DEG 2 TIME 8s

> make && ./main 50 1000 1
clang++ -O3 -march=native main.cpp -o main -I./pml/ntl-extras/ -I./pml/ntl-extras/util/ ./pml/ntl-extras/lib/*.o -lntl -pthread
N 50 DEG 1000 TERMS 1
Using prime 28407454060060787
o: 102002
non pivot 50
[T] SingleApproximant main loop: 3302 ms
[T] SingleApproximant full: 3399 ms
[T] Padé approximation: 1223 ms
[i] Numerator/denominator degrees 50000 / 50000
[i] Result: degree 51000 coeffs 1 6830999029564501 5216029050962576 14445825748146259 ... 28013702577702058
[T] Our total time: 6281 ms

[T] determinant_via_linsolve: 10327 ms
[i] Result: degree 51000: 1 6830999029564501 5216029050962576 14445825748146259 ... 28013702577702058

[T] determinant_via_evaluation_geometric: 21125 ms
[i] Result: degree 51000: 1 6830999029564501 5216029050962576 14445825748146259 ... 28013702577702058

[T] determinant_generic_knowing_degree: 18467 ms
[i] Result: degree 51000: 1 6830999029564501 5216029050962576 14445825748146259 ... 28013702577702058
*/

#include <util/util.h>
#include <mat_lzz_pX_extra/mat_lzz_pX_utils.h>
#include <mat_lzz_pX_extra/mat_lzz_pX_determinant.h>

#include <lzz_pX_CRT.h>
#include <mat_lzz_pX_extra/mat_lzz_pX_approximant.h>
#include <climits>

// #include <flint/nmod_mpoly.h>
// #include <flint/nmod_poly_mat.h>

#define STRLEN 10000000

#define COMPUTE_HARD_COLUMN 1

using namespace std;
using namespace NTL;
using namespace PML;

static zz_p nonzero_rand_zz_p() {
    zz_p a;
    do { a = random_zz_p(); } while (IsZero(a));
    return a;
}
static zz_pX make_sparse_poly(long DEG, long mid_terms = 3) {
    zz_pX f;

    // Required terms
    SetCoeff(f, 0,   nonzero_rand_zz_p());  // degree 0 term
    SetCoeff(f, DEG, nonzero_rand_zz_p());  // degree DEG term (ensures deg(f) == DEG)

    // Pick distinct intermediate degrees
    unordered_set<long> used;
    used.reserve((size_t)mid_terms * 2);
    used.insert(0);
    used.insert(DEG);

    while ((long)used.size() < mid_terms + 2) {
        long d = 1 + RandomBnd(DEG); // 1..DEG-1
        used.insert(d);
    }

    // Assign random nonzero coefficients to the intermediate degrees
    for (long d : used) {
        if (d == 0 || d == DEG) continue;
        SetCoeff(f, d, nonzero_rand_zz_p());
    }

    return f;
}


struct Term {
    long int exp;
    zz_p coef;
    long int idx;

    bool operator<(const Term &other) const {
        return exp < other.exp;
    }
};

static vector<long> pivot_rows_row_echelon(const Mat<zz_p>& A, long want_rank) {
    const long m = A.NumRows();
    const long n = A.NumCols();

    Mat<zz_p> M = A;                         // working copy
    vector<long> row_id(m);
    for (long i = 0; i < m; ++i) row_id[i] = i;

    vector<long> pivots;
    pivots.reserve(want_rank);

    long r = 0;
    for (long c = 0; c < n && r < m; ++c) {
        long piv = r;
        while (piv < m && IsZero(M[piv][c])) ++piv;
        if (piv == m) continue;

        if (piv != r) {
            Vec<zz_p> tmp = M[piv];
            M[piv] = M[r];
            M[r] = tmp;
            swap(row_id[piv], row_id[r]);
        }

        // Normalize pivot row
        zz_p inv_p = inv(M[r][c]);
        for (long j = c; j < n; ++j) M[r][j] *= inv_p;

        // Eliminate below
        for (long i = r + 1; i < m; ++i) {
            if (IsZero(M[i][c])) continue;
            zz_p f = M[i][c];
            for (long j = c; j < n; ++j) M[i][j] -= f * M[r][j];
        }

        pivots.push_back(row_id[r]);
        ++r;

        if ((long)pivots.size() == want_rank) break;
    }
    return pivots;
}
// Mat<zz_p> transpose(const Mat<zz_p>& A) {
//     Mat<zz_p> B;
//     B.SetDims(A.NumCols(), A.NumRows());
//     for (long i = 0; i < A.NumRows(); ++i)
//         for (long j = 0; j < A.NumCols(); ++j)
//             B[j][i] = A[i][j];
//     return B;
// }
// Returns an m-vector of polynomials (zz_pX) like your function.
Vec<zz_pX> SingleApproximantOpt_NTL(const Mat<zz_pX>& G, long o) {
    bool progress = false;
    const long m = G.NumRows();
    const long n = G.NumCols() - 1;

    if (m <= n) throw runtime_error("Need m > n (typically m = n+1).");
    if (o <= 0) throw runtime_error("Order o must be positive.");

    // printf("step 1\n");
    
    // Build constant-term matrix Gc (m x n)
    Mat<zz_p> Gc;
    Gc.SetDims(m, n);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
            Gc[i][j] = coeff(G[i][j], 0);

    // printf("step 2\n");

    // Find pivot rows (need rank n)
    vector<long> pivots = pivot_rows_row_echelon(Gc, n);
    if ((long)pivots.size() != n)
        throw runtime_error("Constant matrix Gc does not have full column rank n.");

    // printf("step 3\n");

    // Build G0 = Gc[pivots, :]
    Mat<zz_p> G0;
    G0.SetDims(n, n);
    for (long r = 0; r < n; ++r) {
        // cout << "pivot " << r << " " << pivots[r] << endl;
        for (long c = 0; c < n; ++c)
            G0[r][c] = Gc[pivots[r]][c];
    }

    // printf("step 4\n");

    // Invert G0
    Mat<zz_p> G0_inv;
    inv(G0_inv, G0);
    Mat<zz_p> G0_inv_T = transpose(G0_inv);
        // throw runtime_error("G0 is not invertible (unexpected if pivots chosen correctly).");

    // printf("step 5\n");

    // Find first non-pivot row i
    vector<char> is_pivot(m, 0);
    for (long pr : pivots) is_pivot[pr] = 1;

    // printf("step 6\n");

    long i_nonpivot = -1;
    for (long i = 0; i < m; ++i) {
        if (!is_pivot[i]) { i_nonpivot = i; break; }
    }
    if (i_nonpivot < 0)
        throw runtime_error("No non-pivot row found (need m>n).");

    cout << "non pivot " << i_nonpivot << endl;
    // printf("step 7\n");
    // row = Gc[i_nonpivot, :]
    Vec<zz_p> row;
    row.SetLength(n);
    for (long j = 0; j < n; ++j) row[j] = Gc[i_nonpivot][j];

    // printf("step 8\n");
    // z0 = - row * G0_inv   (1xn)
    Vec<zz_p> z0;
    z0.SetLength(n);
    for (long j = 0; j < n; ++j) {
        zz_p s = zz_p(0);
        // for (long k = 0; k < n; ++k) s += row[k] * G0_inv[k][j];
        for (long k = 0; k < n; ++k) s += row[k] * G0_inv_T[j][k];
        z0[j] = -s;
    }

    // printf("step 9\n");
    // inds = pivots + [i_nonpivot]
    vector<long> inds = pivots;
    inds.push_back(i_nonpivot);

    // Pre-extract sparse terms for each polynomial G[i][j], restricted to exponents < o
    vector<vector<vector<Term>>> Gsparse(m, vector<vector<Term>>(n));
    vector<vector<Term>> GsparseT(n);

    for (long r = 0; r < m; ++r) {
        for (long c = 0; c < n; ++c) {
            long dmax = min(o - 1, deg(G[r][c]));
            auto& terms = Gsparse[r][c];
            terms.reserve(8);
            for (long e = 0; e <= dmax; ++e) {
                zz_p a = coeff(G[r][c], e);
                if (!IsZero(a)) {
                    terms.push_back({e, a, r});
                    GsparseT[c].push_back(terms.back());
                }
            }
        }
    }
    for (long c = 0; c < n; ++c) {
        sort(GsparseT[c].begin(), GsparseT[c].end());
    }

    // zcoeff: (n+1) polynomials, each stored as coefficient array length o
    // zcoeff[0..n-1][0] = z0, zcoeff[n][0] = 1
    vector<Vec<zz_p>> zcoeff(n + 1);
    for (long k = 0; k < n + 1; ++k) {
        zcoeff[k].SetLength(o);
        for (long t = 0; t < o; ++t) zcoeff[k][t] = zz_p(0);
    }
    for (long k = 0; k < n; ++k) zcoeff[k][0] = z0[k];
    zcoeff[n][0] = zz_p(1);

    
    Vec<zz_p> v;
    v.SetLength(n);
    Vec<zz_p> t;
    t.SetLength(n);

    // cout << "terms  " << Gsparse[0][0].size() << endl;
    // cout << "termsT " << GsparseT[0].size() << endl;

    auto t0 = chrono::steady_clock::now();
    // Main loop: d = 1..o-1
    for (long d = 1; d < o; ++d) {
        if (progress && (d % 100 == 0))
            cerr << "d=" << d << "/" << (o - 1) << "\n";

        // v = coefficient of x^d in (z(x) * G(x)) for each column
        for (long j = 0; j < n; ++j) {
            zz_p acc = zz_p(0);
            for (const auto& tm : GsparseT[j]) {
                if (tm.exp <= d) {
                    acc += tm.coef * zcoeff[tm.idx][d - tm.exp];
                }
                else {
                    break;
                }
            }
            v[j] = acc;
        }

        // t = v * G0_inv  (1xn)
        for (long j = 0; j < n; ++j) {
            zz_p s = zz_p(0);
            // for (long k = 0; k < n; ++k) s += v[k] * G0_inv[k][j];
            for (long k = 0; k < n; ++k) {
                s += v[k] * G0_inv_T[j][k];
            }
            zcoeff[j][d] -= s;
        }
    }
    auto t1 = chrono::steady_clock::now();
    cerr << "[T] SingleApproximant main loop: " << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() << " ms\n";

    // Build result vector of polynomials length m (only inds positions are nonzero)
    Vec<zz_pX> res;
    res.SetLength(m);
    // for (long r = 0; r < m; ++r) clear(res[r]);

    for (long zi = 0; zi < n+1; ++zi) {
        long pi = inds[zi];
        zz_pX f;
        // clear(f);
        for (long e = 0; e < o; ++e) {
            if (!IsZero(zcoeff[zi][e])) SetCoeff(f, e, zcoeff[zi][e]);
        }
        res[pi] = f;
    }

    return res;
}

struct Mat2X {
    zz_pX a00, a01, a10, a11; // [[a00 a01],[a10 a11]]

    static Mat2X Identity() {
        Mat2X M;
        clear(M.a00); clear(M.a01); clear(M.a10); clear(M.a11);
        SetCoeff(M.a00, 0, 1);
        SetCoeff(M.a11, 0, 1);
        return M;
    }

    static Mat2X Swap() { // [[0 1],[1 0]]
        Mat2X M;
        clear(M.a00); clear(M.a11);
        clear(M.a01); clear(M.a10);
        SetCoeff(M.a01, 0, 1);
        SetCoeff(M.a10, 0, 1);
        return M;
    }

    // Matrix multiply: this * B
    Mat2X operator*(const Mat2X& B) const {
        Mat2X C;
        C.a00 = a00*B.a00 + a01*B.a10;
        C.a01 = a00*B.a01 + a01*B.a11;
        C.a10 = a10*B.a00 + a11*B.a10;
        C.a11 = a10*B.a01 + a11*B.a11;
        return C;
    }

    // Apply to (x,y): returns (u,v) = (a00*x+a01*y, a10*x+a11*y)
    void apply(zz_pX& u, zz_pX& v, const zz_pX& x, const zz_pX& y) const {
        u = a00*x + a01*y;
        v = a10*x + a11*y;
    }
};

// Truncate to degree < n  (i.e., mod x^n)
static zz_pX trunc_mod_xn(const zz_pX& f, long n) {
    zz_pX g;
    clear(g);
    long d = deg(f);
    if (d < 0) return g;
    long top = min(d, n-1);
    for (long i = 0; i <= top; ++i) {
        zz_p c = coeff(f, i);
        if (!IsZero(c)) SetCoeff(g, i, c);
    }
    return g;
}

static void HalfGCD_for_Pade(const zz_pX& f_in, long n, zz_pX* num, zz_pX* den) {
    zz_pX xn;
    clear(xn);
    SetCoeff(xn, n, 1);

    zz_pX f = trunc_mod_xn(f_in, n);
    zz_pXMatrix M;
    long d_red = n/2 + 1;
    HalfGCD(M, xn, f, d_red);

    zz_pX r0 = M(0,0) * xn + M(0,1) * f;
    zz_pX r1 = M(1,0) * xn + M(1,1) * f;
    if (num) *num = r1;
    if (den) *den = M(1,1);
}

Vec<zz_pX> mul_scalar(const Vec<zz_pX>& v, const zz_pX& s) {
    Vec<zz_pX> out;
    out.SetLength(v.length());
    for (long i = 0; i < v.length(); ++i)
        out[i] = v[i] * s;
    return out; // RVO / move
}
zz_pX dot(const Vec<zz_pX>& a, const Vec<zz_pX>& b) {
    if (a.length() != b.length())
        throw runtime_error("dot: length mismatch");

    zz_pX s;
    clear(s);
    for (long i = 0; i < a.length(); ++i)
        s += a[i] * b[i];

    return s; // polynomial
}
Vec<zz_pX> get_col(const Mat<zz_pX>& M, long j) {
    if (j < 0 || j >= M.NumCols())
        throw runtime_error("get_col: column index out of range");

    Vec<zz_pX> col;
    col.SetLength(M.NumRows());
    for (long i = 0; i < M.NumRows(); ++i)
        col[i] = M[i][j];

    return col;
}

bool CheckSingleApproximant(const Mat<zz_pX>& G, const Vec<zz_pX>& z, long o, bool check_has_unit=true) {
    const long m = G.NumRows();
    const long n = G.NumCols() - 1;
    if (z.length() != m) throw runtime_error("z length != #rows(G)");

    // Check z*G ≡ 0 (mod x^o)
    for (long j = 0; j < n; ++j) {
        zz_pX acc; clear(acc);
        for (long i = 0; i < m; ++i) acc += z[i] * G[i][j];

        for (long e = 0; e < o; ++e) {
            if (!IsZero(coeff(acc, e))) {
                cerr << "FAIL: column " << j << " has nonzero coeff at x^" << e << "\n";
                return false;
            }
        }
    }

    // Optional: check some entry equals the constant polynomial 1
    if (check_has_unit) {
        bool found = false;
        for (long i = 0; i < m; ++i) {
            if (deg(z[i]) == 0 && coeff(z[i], 0) == zz_p(1)) { found = true; break; }
        }
        if (!found) {
            cerr << "WARN: no coordinate is exactly the constant polynomial 1\n";
            // not necessarily fatal if you changed the construction
        }
    }

    return true;
}

struct SparseTerm {
    zz_p ratio;
    zz_p cur;
};

static void determinant_via_evaluation_geometric_sparse(zz_pX & det, const Mat<zz_pX> & pmat) {
    const long d = deg(pmat);
    const long dim = pmat.NumRows();
    const long nb_points = dim * d + 1;

    zz_pX_Multipoint_Geometric ev = get_geometric_points(nb_points);
    const zz_p q = ev.get_q();
    const zz_p s = ev.get_s();

    // Precompute sparse terms for each entry with geometric progression across points.
    auto t_pre0 = chrono::steady_clock::now();
    vector<vector<SparseTerm>> terms_flat(dim * dim);
    for (long i = 0; i < dim; ++i) {
        for (long j = 0; j < dim; ++j) {
            const zz_pX& f = pmat[i][j];
            long df = deg(f);
            if (df < 0) continue;
            auto &tv = terms_flat[i * dim + j];
            tv.reserve(8);
            for (long e = 0; e <= df; ++e) {
                zz_p a = coeff(f, e);
                if (IsZero(a)) continue;
                zz_p qs = power(q, e);
                zz_p ss = power(s, e);
                zz_p cs = a * ss;
                tv.push_back({qs, cs}); // ratio = q^e, cur = coeff * s^e
            }
        }
    }
    auto t_pre1 = chrono::steady_clock::now();
    cerr << "[T] sparse_eval precompute: "
         << chrono::duration_cast<chrono::milliseconds>(t_pre1 - t_pre0).count()
         << " ms\n";

    Vec<zz_p> det_evals(INIT_SIZE, nb_points);
    Mat<zz_p> eval;
    eval.SetDims(dim, dim);

    long sz = dim * dim;
    for (long k = 0; k < nb_points; ++k) {
        long ii = 0, jj = 0;
        for (long i = 0; i < sz; ++i) {
            zz_p sum = zz_p(0);
            auto &tv = terms_flat[i];
            for (auto &t : tv) {
                sum += t.cur;
            }
            eval[ii][jj] = sum;
            jj++;
            if (jj >= dim) {
                jj -= dim;
                ii++;
            }

        }

        determinant(det_evals[k], eval);

        // Advance all term currents for next point: cur *= ratio
        for (long i = 0; i < sz; ++i) {
            auto &tv = terms_flat[i];
            for (auto &t : tv) {
                t.cur *= t.ratio;
            }
        }
    }

    ev.interpolate(det, det_evals);
}

int main(int argc, char **argv)
{
    if (argc != 4) {
        printf("Usage: %s N DEG TERMS\n", argv[0]);
        return 0;
    }
    const long N   = atoi(argv[1]);
    const long DEG = atoi(argv[2]);
    const long TERMS = atoi(argv[3]);

    cout << "N " << N << " DEG " << DEG << " TERMS " << TERMS << endl;

    const long long p = 28407454060060787ll; 
    zz_p::init(p);

    cout << "Using prime " << p << endl;

    // Seed NTL RNG (optional)
    SetSeed(ZZ((long)12345));

    Mat<zz_pX> A;
    A.SetDims(N+1, N+1);

    for (long i = 0; i < N+1; ++i) {
        for (long j = 0; j < N+1; ++j) {
            A[i][j] = make_sparse_poly(DEG, TERMS);
            // cout << A[i][j][0] << " ";
        }
        // cout << endl;
    }

    // Example: print one entry (remove if you don't want output)
    // cout << "A[0][0] = " << A[0][0] << "\n";
    // cout << "deg(A[0][0]) = " << deg(A[0][0]) << "\n";

    long o = DEG * (N + 1) * 2 + 2;

    cout << "o: " << o << endl;

    auto t00 = chrono::steady_clock::now();
    auto t0 = chrono::steady_clock::now();
    Vec<zz_pX> vec = SingleApproximantOpt_NTL(A, o);
    auto t1 = chrono::steady_clock::now();
    cerr << "[T] SingleApproximant full: " << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() << " ms\n";

    // cout << "approximant ok? " << CheckSingleApproximant(A, vec, o) << endl;
    // cout << deg(vec[0]) << " " << deg(vec[1]) << endl << endl;
    // cout << "my res: " << res << endl;

    zz_pX num, den;
    t0 = chrono::steady_clock::now();
    HalfGCD_for_Pade(vec[0], o, &num, &den);
    t1 = chrono::steady_clock::now();
    cerr << "[T] Padé approximation: " << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() << " ms\n";

    cout << "[i] Numerator/denominator degrees " << deg(num) << " / " << deg(den) << endl;
    // cout << "EXPECT " << num << endl;
    // cout << "   GOT " << trunc_mod_xn(den * vec[0], o) << endl << endl;

    // cout << "num" << num << endl;
    // for (int i = 0; i < N+1; i++) {
        // cout << "vec? " << deg(vec[i]) << endl;// << " " <<  ker[0] << endl << endl;
    // }

    Vec<zz_pX> ker = mul_scalar(vec, den);
    for (int i = 0; i < ker.length(); i++) {
        ker[i] = trunc_mod_xn(ker[i], o);
    }
    // cout << "ker?" << ker << endl << endl;
    // for (int i = 0; i < N+1; i++) {
        // cout << "ker? " << deg(ker[i]) << endl;// << " " <<  ker[0] << endl << endl;
    // }
    // cout << endl;
    
    // cout << "ker * A degs: ";
    // for (int i = 0; i < N ; i++) {
    //     cout << deg(dot(get_col(A, i), ker)) << " ";
    // }
    // cout << endl;

    zz_pX res = dot(get_col(A, N), ker);
    // cout << "ker? " << deg(ker[1]) << endl << endl;// << " " << ker[1] << endl << endl;
    auto ires0 = inv(res[0]);
    for (int i = deg(res); i >= 0; i--) {
        res[i] *= ires0;
    }
    // cout << endl;
    cout << "[i] Result: degree " << deg(res) << " coeffs " << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << " ... " << res[deg(res)] << endl;
    t1 = chrono::steady_clock::now();
    uint64_t our_elapsed = chrono::duration_cast<chrono::milliseconds>(t1 - t00).count();
    cerr << "[T] Our total time: " << our_elapsed << " ms\n";
    cout << endl;
    {
        zz_pX det;

        if (1) {
            t0 = chrono::steady_clock::now();
            determinant_via_linsolve(det, A);
            t1 = chrono::steady_clock::now();
            cerr << "[T] determinant_via_linsolve: " << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() << " ms\n";
            for (int i = deg(det); i >= 0; i--) {
                det[i] /= det[0];
            }
            cout << "[i] Result: degree " << deg(det) << ": " << det[0] << " " << det[1] << " " << det[2] << " " << det[3] << " ... " << det[deg(det)] << endl;
            cout << endl;
        }

        if (1) {
            t0 = chrono::steady_clock::now();
            determinant_via_evaluation_geometric_sparse(det, A);
            t1 = chrono::steady_clock::now();
            cerr << "[T] determinant_via_evaluation_geometric_sparse: " << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() << " ms\n";
            for (int i = deg(det); i >= 0; i--) {
                det[i] /= det[0];
            }
            cout << "[i] Result: degree " << deg(det) << ": " << det[0] << " " << det[1] << " " << det[2] << " " << det[3] << " ... " << det[deg(det)] << endl;
            cout << endl;
        }

        if (1) {
            t0 = chrono::steady_clock::now();
            determinant_via_evaluation_geometric(det, A);
            t1 = chrono::steady_clock::now();
            cerr << "[T] determinant_via_evaluation_geometric: " << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() << " ms\n";
            for (int i = deg(det); i >= 0; i--) {
                det[i] /= det[0];
            }
            cout << "[i] Result: degree " << deg(det) << ": " << det[0] << " " << det[1] << " " << det[2] << " " << det[3] << " ... " << det[deg(det)] << endl;
            cout << endl;
        }

        if (0) {
            t0 = chrono::steady_clock::now();
            determinant_via_evaluation_FFT(det, A);
            t1 = chrono::steady_clock::now();
            cerr << "[T] determinant_via_evaluation_FFT: " << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() << " ms\n";
            for (int i = deg(det); i >= 0; i--) {
                det[i] /= det[0];
            }
            cout << "[i] Result: degree " << deg(det) << ": " << det[0] << " " << det[1] << " " << det[2] << " " << det[3] << " ... " << det[deg(det)] << endl;
            cout << endl;
        }

        if (1) {
            t0 = chrono::steady_clock::now();
            determinant_generic_knowing_degree(det, A, DEG * (N+1));
            t1 = chrono::steady_clock::now();
            cerr << "[T] determinant_generic_knowing_degree: " << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() << " ms\n";
            for (int i = deg(det); i >= 0; i--) {
                det[i] /= det[0];
            }
            cout << "[i] Result: degree " << deg(det) << ": " << det[0] << " " << det[1] << " " << det[2] << " " << det[3] << " ... " << det[deg(det)] << endl;
            cout << endl;
        }

        // cout << "coef0 " << det[0] << endl;
        
        // cout << res << endl;
        // cout << det << endl;
    }   

    return 0;
}