# Magic Pot: Cryptanalysis of Full AIM2 in the Standard and Related-/reused-Key Settings Using New Elimination Framework

This repository contains supporting code for the [EUROCRYPT 2026](https://eurocrypt.iacr.org/2026/) paper ([doi](10.1007/978-3-032-25333-0_3)], ([eprint full version](https://eprint.iacr.org/2026/xx)):

> **Magic Pot: Cryptanalysis of Full AIM2 in the Standard and Related-/reused-Key Settings Using New Elimination Framework**

by Alex Biryukov, Pablo García Fernández & Aleksei Udovenko.

The work was funded by Luxembourg's FNR projects CryptoFin (C22/IS/17415825), PQseal (C24/IS/18978392), NCER22/IS/16570468/NCER-FT, and FinnovationHub funded by the Ministry of Finance, Government of Luxembourg.

```bib
@InProceedings{10.1007/978-3-032-25333-0_3,
    author="Biryukov, Alex and Garc{\'i}a Fern{\'a}ndez, Pablo and Udovenko, Aleksei",
    editor="Daemen, Joan and Thom{\'e}, Emmanuel",
    title="Magic Pot: Cryptanalysis of Full AIM2 in the Standard and Related-/reused-Key Settings Using New Elimination Framework",
    booktitle="Advances in Cryptology -- EUROCRYPT 2026",
    year="2026",
    publisher="Springer Nature Switzerland",
    address="Cham",
    pages="63--90",
    isbn="978-3-032-25333-0",
    doi="10.1007/978-3-032-25333-0_3"
}
```

A copy of this repository is available at [zenodo.org](https://doi.org/10.5281/zenodo.18727305).



## Setup
Requirements: Rust setup, maturin, SageMath.

Installation of the ff-rust library:

```sh
$ sage -python -m pip install ./ff_rust
```

Or in a notebook:

```py
%pip install ./ff_rust
```

## Files

- [aim2.py](./aim2.py) contains SageMath implementation of AIM2, compatible with the reference implementations.
- [attack-aim2.py](./attack-aim2.py) contains SageMath implementation of the estimation and key recovery attack.
- [polynomial.py](./polynomial.py) contains generic (unoptimized) multivariate polynomial ring implementation, to track the equations.
- [halfgcd.py](./halfgcd.py) contains univariate polynomial arithmetic wrappers (gcd, halfgcd, division/remainder) to use FAFFT for fast operations.
- [keypairs_from_reference_implementation](./keypairs_from_reference_implementation) folder contains C program to generate (iv,ct) pairs from a reference C implementation, as well as pre-generated sample files. These files contain the secret key for verification purposes, but the key can be erased to make sure the key recovery is not cheating.

## Attack

```sh
$ sage attack-aim2.py -h
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
```

The parameters should be clear from the usage help. There are two main options of obtaining (iv,ct) pairs:

**Option 1:** generate random IVs to construct the system. ell is determined automatically, L=1 by default.

```sh
$ sage attack-aim2.py -n 128

$ sage attack-aim2.py -n 128 -L 50  # few seconds

$ sage attack-aim2.py -n 128 -L 15  # 1 hour

$ sage attack-aim2.py -n 128 -l 5  # override with 5 branches
```

**Option 2:** run on file with (iv, ct) pairs generated from C reference implementation (see [simple_keypairs.c](keypairs_from_reference_implementation/simple_keypairs.c)). Parameters are determined automatically, number of IVs can be overriden to a smaller quantity.

```sh
$ sage attack-aim2.py --pairs-file keypairs_from_reference_implementation/keypairs128.txt -L 30  # out of 100 in the file

$ sage attack-aim2.py --pairs-file keypairs_from_reference_implementation/keypairs256.txt
```


## Example

```sh
$ sage attack-aim2.py -n 128 -L 30 | tee logs/log_n128_ell2_L30.txt
[i] Seed: 3176760365509728981
[i] Secret key: GF128Element(27d6539e87ba3458344d07ad8988b912)

Generating (iv, ct) pairs...
100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 30/30 [00:00<00:00, 248.46it/s]


[i] AIM/attack parameters:
- n = 128 bits
- ell = 2 branches
- L = 30 IVs
- M = 8 squarings (1 ... 2^(M-1) 2^M)
- W = 257 t-monomials
- Only estimate ? False

[I] Time complexity rough estimate = 2^38.39 field ops

Generating base equations... (using linearized polynomials)
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 30/30 [00:02<00:00, 11.42it/s]
0it [00:00, ?it/s]
[i] Initial equations: 32
[i] Squared equations: 288
[i] W = t-monomials in equations: 257

Arranging equations...
288it [00:01, 200.17it/s]
single_approximant: len(flat_coeffs) = 122931
single_approximant: d = 1683 / 16834 (10.0 %)
single_approximant: d = 3366 / 16834 (20.0 %)
single_approximant: d = 5049 / 16834 (30.0 %)
single_approximant: d = 6732 / 16834 (40.0 %)
single_approximant: d = 8415 / 16834 (50.0 %)
single_approximant: d = 10098 / 16834 (60.0 %)
single_approximant: d = 11781 / 16834 (70.0 %)
single_approximant: d = 13464 / 16834 (80.0 %)
single_approximant: d = 15147 / 16834 (90.0 %)
single_approximant: d = 16830 / 16834 (100.0 %)

[i] Full equation matrix F: H=288 x W=257

Computing pivots (p=0)
[i] Rank = #pivots = 257/257   fingerprint 2cb06969a1aa15e7
[I] Row-degree sum ξ: experimental = 2^13.04, estimated = 2^13.84
[I] Total nonzero coefficients 123411 <= 2W^2 = 132098 ? True
[I] Time complexity 2^31.90 field ops = 2^24.90 encs  (1/n convention)
[I] Memory complexity : 2^22.25 field elems = 2^26.25 bytes = 0.07 GB
[R1] n=128  ell= 2  L= 30  M=  8  W= 257  ξ=2^ 13.04  Time=2^ 31.90 field ops  Mem=2^ 22.25 field elems
[R2] n=128  ell= 2  L= 30  M=  8  W= 257  ξ=2^ 13.04  Time=2^ 24.90 enc  Mem=2^ 26.25 bytes

[i] Running full attack 0.1 GB needed, 26.5 GB available

[o] Approximant order o = 16834 = 2^14.04
[i] FAFFT max degree 2^16
[T] Main attack part, start time measurement  | memory usage: 0.26 GB

[.] Computing pivots of G0 for inversion...
[.] Inversion of base matrix...
[.] Single approximant...
[T] TIME after single approximant: 30.03308129310608  | memory usage: 0.33 GB

[i] Kernel len = 257 degrees = [16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 16833, 0]
[T] TIME after first HalfGCD: 30.83499312400818  | memory usage: 0.33 GB
[i] Degrees of kernel * F [-1, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8415, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8414, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8412, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8408, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8400, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8384, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8352, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8288, 8160]
[i] Got univariate poly in p, degree 8416
[T] TIME to get univariate polynomial: 34.36897897720337  | memory usage: 0.27 GB

[T] TIME pre gcd (x^(2^n)-x mod f): 39.83170938491821  | memory usage: 0.27 GB
gcd polynomial 203496734187870837048844359436221223931*x^3 + 175969891944308531575527162406742962215*x^2 + 67431397192295501745274199971634941019*x + 237055397823692163691616443013404842587
[T] TIME post gcd: 40.66392731666565  | memory usage: 0.27 GB
found roots:
- 1428750bca12e1b2063ab4527764dbcf | 0 / 30 match
- 27d6539e87ba3458344d07ad8988b912 | 30 / 30 match
- 80dfc4b2e7c90a6040c8c7222f4bb8bc | 0 / 30 match
[T] TIME full attack (main part only): 41.01423716545105  | memory usage: 0.27 GB
[i] Real secret key for comparison: GF128Element(27d6539e87ba3458344d07ad8988b912)
```

## Reproducing Table 4

[Table4.sh](./Table4.sh) runs estimations and attacks for Table 4 (and a bit more) in a few hours on a laptop, where most of the time is spent on the three practical key recovery attacks with small number of IVs. Using grep we can easily obtain the summary. Lines with `[R1]` contain summary with field operations & elements units, lines with `[R2]` contain summary with encryption & bytes units. Full logs are available in the [./logs](./logs) folder and contain extra information, such as running time, real memory consumption, etc.

```
$ bash Table4.sh
...
$ grep '\[R1]' logs/n*.log -h | sort | tee Table4.log
[R1] n=128  ell= 2  L=127  M=  1  W= 257  ξ=2^  8.59  Time=2^ 27.30 field ops  Mem=2^ 19.71 field elems
[R1] n=128  ell= 2  L= 15  M= 15  W= 257  ξ=2^ 19.25  Time=2^ 38.44 field ops  Mem=2^ 28.26 field elems
[R1] n=128  ell= 2  L=  2  M= 64  W= 257  ξ=2^ 66.32  Time=2^ 88.64 field ops  Mem=2^ 75.33 field elems
[R1] n=128  ell= 2  L=  1  M= 85  W= 257  ξ=2^ 87.32  Time=2^110.66 field ops  Mem=2^ 96.33 field elems
[R1] n=128  ell= 3  L=  1  M= 96  W= 385  ξ=2^ 98.32  Time=2^122.26 field ops  Mem=2^107.91 field elems
[R1] n=128  ell= 4  L=  1  M=102  W= 513  ξ=2^105.00  Time=2^129.32 field ops  Mem=2^115.00 field elems
[R1] n=128  ell= 4  L=  1  M=102  W= 513  ξ=2^105.00  Time=2^129.32 field ops  Mem=2^115.00 field elems
[R1] n=128  ell= 5  L=  1  M=106  W= 641  ξ=2^109.46  Time=2^134.06 field ops  Mem=2^119.78 field elems
[R1] n=128  ell= 6  L=  1  M=109  W= 769  ξ=2^112.70  Time=2^137.52 field ops  Mem=2^123.29 field elems
[R1] n=128  ell= 7  L=  1  M=112  W= 897  ξ=2^115.17  Time=2^140.18 field ops  Mem=2^125.98 field elems
[R1] n=128  ell= 8  L=  1  M=113  W=1025  ξ=2^117.09  Time=2^142.26 field ops  Mem=2^128.09 field elems
[R1] n=128  ell= 9  L=  1  M=115  W=1153  ξ=2^118.70  Time=2^144.03 field ops  Mem=2^129.87 field elems

[R1] n=192  ell= 2  L=191  M=  1  W= 385  ξ=2^  9.17  Time=2^ 29.02 field ops  Mem=2^ 20.88 field elems
[R1] n=192  ell= 2  L= 31  M= 11  W= 385  ξ=2^ 16.78  Time=2^ 36.83 field ops  Mem=2^ 26.39 field elems
[R1] n=192  ell= 2  L=  2  M= 96  W= 385  ξ=2^ 98.32  Time=2^122.32 field ops  Mem=2^107.91 field elems
[R1] n=192  ell= 2  L=  1  M=128  W= 385  ξ=2^130.00  Time=2^155.04 field ops  Mem=2^139.59 field elems
[R1] n=192  ell= 3  L=  1  M=144  W= 577  ξ=2^146.32  Time=2^171.97 field ops  Mem=2^156.49 field elems
[R1] n=192  ell= 4  L=  1  M=153  W= 769  ξ=2^156.17  Time=2^182.19 field ops  Mem=2^166.76 field elems
[R1] n=192  ell= 5  L=  1  M=160  W= 961  ξ=2^162.81  Time=2^189.10 field ops  Mem=2^173.72 field elems
[R1] n=192  ell= 6  L=  1  M=164  W=1153  ξ=2^167.58  Time=2^194.09 field ops  Mem=2^178.76 field elems
[R1] n=192  ell= 7  L=  1  M=168  W=1345  ξ=2^171.17  Time=2^197.85 field ops  Mem=2^182.56 field elems
[R1] n=192  ell= 8  L=  1  M=170  W=1537  ξ=2^174.00  Time=2^200.84 field ops  Mem=2^185.59 field elems
[R1] n=192  ell= 9  L=  1  M=172  W=1729  ξ=2^176.25  Time=2^203.23 field ops  Mem=2^188.00 field elems

[R1] n=256  ell= 3  L=382  M=  1  W= 769  ξ=2^ 10.17  Time=2^ 31.97 field ops  Mem=2^ 22.87 field elems
[R1] n=256  ell= 3  L= 94  M=  7  W= 769  ξ=2^ 14.54  Time=2^ 36.38 field ops  Mem=2^ 25.34 field elems
[R1] n=256  ell= 3  L=  2  M=153  W= 769  ξ=2^156.17  Time=2^182.23 field ops  Mem=2^166.76 field elems
[R1] n=256  ell= 3  L=  1  M=192  W= 769  ξ=2^194.32  Time=2^221.19 field ops  Mem=2^204.91 field elems
[R1] n=256  ell= 4  L=  1  M=204  W=1025  ξ=2^207.32  Time=2^234.55 field ops  Mem=2^218.32 field elems
[R1] n=256  ell= 5  L=  1  M=213  W=1281  ξ=2^216.17  Time=2^243.67 field ops  Mem=2^227.49 field elems
[R1] n=256  ell= 6  L=  1  M=219  W=1537  ξ=2^222.46  Time=2^250.17 field ops  Mem=2^234.05 field elems
[R1] n=256  ell= 7  L=  1  M=224  W=1793  ξ=2^227.17  Time=2^255.06 field ops  Mem=2^238.98 field elems
[R1] n=256  ell= 8  L=  1  M=227  W=2049  ξ=2^230.91  Time=2^258.95 field ops  Mem=2^242.91 field elems
[R1] n=256  ell= 9  L=  1  M=230  W=2305  ξ=2^233.91  Time=2^262.08 field ops  Mem=2^246.08 field elems
```