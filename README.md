# AK_Riemann_zeta_gamma_C

## Overview

This is a C language integration and demo of the Riemann_zeta and ln_gamma
functions that Dr. Alexey Kuznetsov provided on
[GitHub](https://github.com/Alexey-Kuznetsov-math/Riemann_zeta-and-Gamma),
specifically, the quad precision Fortran versions.  They are usable in C++ as
extern "C" functions.

A GNU environment (GNU make, gcc, gfortran, GNU libc, etc.) is assumed.  The
source code of included programs uses GNU extensions.

My additions are © 2026 David Flater but are provided under the same terms
as the originals (BSD 3-Clause License).

## Functional profile

You may read Dr. Kuznetsov's own overview
[here](https://github.com/Alexey-Kuznetsov-math/Riemann_zeta-and-Gamma).

Background:  [Riemann
zeta](https://en.wikipedia.org/wiki/Riemann_zeta_function) ζ(s) and
[gamma](https://en.wikipedia.org/wiki/Gamma_function) Γ(z) are described on
Wikipedia.  [[Gourdon and Sebah
2003](http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.pdf)]
provides an introduction to previous algorithms for approximating zeta in
regions where a Dirichlet series summation does not suffice (Euler–Maclaurin,
Borwein, Riemann–Siegel, Odlyzko–Schönhage).  See also [[Borwein
1995](https://www.cecm.sfu.ca/~pborwein/PAPERS/P155.pdf)] for Borwein
Algorithms 2 and 3, which Gourdon and Sebah referred to as Proposition 1 and
2 under the convergence of alternating series method.

The function Riemann_zeta implements a new algorithm, hereinafter called
Kuznetsov 2025, that is described in
[[preprint](https://arxiv.org/abs/2503.09519)].  It combines the main sum of
Riemann–Siegel with a simple approximation of the remainder term.

Riemann_zeta divides the domain as follows for Re(s) ≥ ½ and uses reflection
for Re(s) < ½:

For |Im(s)| < 200:

- If Re(s) < 4, use Euler–Maclaurin.
- Otherwise, use whichever of Euler–Maclaurin or a Dirichlet series has lower
  computational complexity for s.

For |Im(s)| ≥ 200:

- If Re(s) < 4, use Kuznetsov 2025.
- Otherwise, use whichever of Kuznetsov 2025 or a Dirichlet series has lower
  computational complexity for s.

The expected accuracy is about 31 decimal digits for |Im(s)| < 100, 30
decimal digits for |Im(s)| < 1000, etc., losing one digit with each factor of
10 increase in Im(s).  \[This needs some qualification, as it does not hold
for the trivial zeros as Re(s) becomes a sufficiently large, even negative
integer.\]

The function ln_gamma implements a new alternative to [Lanczos
approximation](https://en.wikipedia.org/wiki/Lanczos_approximation) that is
described in [[preprint](https://arxiv.org/abs/2109.12061), [paywalled
DOI](https://doi.org/10.1016/j.cam.2022.114270)] and
[[preprint](https://arxiv.org/pdf/2508.19095)].  It is expected to maintain
accuracy on the order of 5×10⁻³².

Pragmatically, these functions fill a gap in the math library ecosystem for C/C++:

Library | Version | Quad precision | Complex zeta | Complex gamma
--- | --- | --- | --- | ---
[C standard library](https://en.cppreference.com/w/c/numeric/math.html) | C99 | Yes (as extension) | No | No
[C++ standard library](https://cppreference.net/cpp/numeric.html) | C++17 | Yes (as extension) | No | No
[Boost Math](https://www.boost.org/library/latest/math/) | 1.90.0 | Yes | No | No
[GNU Scientific Library](https://www.gnu.org/software/gsl/) | 2.8 | No | No | Yes

## Code

The functions themselves are provided by [Riemann_zeta.f90](Riemann_zeta.f90)
and [Riemann_zeta.h](Riemann_zeta.h).  You build the Fortran as Fortran and
link it with C/C++ code as the demo programs demonstrate.

Riemann_zeta.f90 is a lightly modified version of
[Fortran/Riemann_zeta.f90](https://github.com/Alexey-Kuznetsov-math/Riemann_zeta-and-Gamma/blob/main/Fortran/Riemann_zeta.f90)
as of [commit
3f1e7b662be8a9302f5a55c7cfa16f65fe557566](https://github.com/Alexey-Kuznetsov-math/Riemann_zeta-and-Gamma/commit/3f1e7b662be8a9302f5a55c7cfa16f65fe557566),
Date:  Fri Dec 5 08:18:03 2025 -0500.

The changes that I made are:

- C language integration:  Converted Riemann_zeta.f90 to module form,
following
[MATLAB_Fortran_mex/Riemann_zeta_module.f90](https://github.com/Alexey-Kuznetsov-math/Riemann_zeta-and-Gamma/blob/main/MATLAB_Fortran_mex/Riemann_zeta_module.f90),
to fix forward function references and missing/duplicated definitions for the
constants, and changed the declarations of Riemann_zeta and ln_gamma to use C
bindings.
- Fix misleading identifier:  Renamed N_12 to N_30 and changed the comment to
  reference zeta_30.
- Fix warnings:  Untabified; deleted unused variables.
- Fix NaN result at origin:  Added a special case for ζ(s) = −½ when |s| <
  10⁻¹⁷.
- Clean up:  Stripped trailing whitespace.

## Demo programs

- zeta:  Test Riemann_zeta from the command line.
Usage: zeta sreal simag
- gamma:  Test ln_gamma from the command line.
Usage: gamma zreal zimag
- test:  Run noninteractive test suite for Riemann_zeta and ln_gamma.  This
is
[test.f90](https://github.com/Alexey-Kuznetsov-math/Riemann_zeta-and-Gamma/blob/main/Fortran/test.f90)
manually converted to C to establish that the results are the same as they
were in Fortran.

To build all three programs normally, just type `make`.

To build all three programs in C++ mode, type `make cxx`.  The C++ binaries
will be called zeta_cxx, gamma_cxx, and test_cxx.

## Sample output

- [test_out_f90.txt](test_out_f90.txt):  Sample output of the original
test.f90 built according to the author's instructions.
- [test_out_c.txt](test_out_c.txt):  Sample output of test.c.

Enough digits are printed to show that the values are the same down to the
last bit.

The max error of Test 3 is expected to change from run to run because it uses
random data.  The reported error is typically on the order of 5×10⁻³², but if
one of the randomly chosen test points happens to be near one of the poles of
the gamma function (zero and the negative integers), it can become infinite.

## Other free implementations of zeta for complex s

- mpmath 1.4.0
[mpmath.zeta](https://www.mpmath.org/doc/current/functions/zeta.html)
(Python):  "The implementation uses three algorithms:  the Borwein algorithm
for the Riemann zeta function when s is close to the real line; the
Riemann-Siegel formula for the Riemann zeta function when s is large
imaginary, and Euler-Maclaurin summation in all other cases."  The Borwein
code in mpmath/functions/zeta.py appears to be commented out.  There's
another implementation in mpmath/libmp/gammazeta.py that's using Borwein
Algorithm 2.
- SciPy 1.17.0
[scipy.special.zeta](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.zeta.html)
(Python) ... subprojects/xsf/include/xsf/zeta.h is C++, Borwein Algorithm 2
plus Euler-Maclaurin ...  "Complex riemann-zeta function implementation based
on [Python implementation written by Matt
Haberland](https://colab.research.google.com/drive/1zMDSAJlXCLRqMMtJ0e9nDGjQ8iZCnmn5?usp=sharing)"
... [unfinished work?](https://github.com/scipy/scipy/issues/14073)
- [Search GitHub](https://github.com/search?q=riemann%20zeta&type=repositories)
