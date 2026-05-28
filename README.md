# AK_Riemann_zeta_gamma_C

## Overview

This is a C language integration and demo of some functions that Dr. Alexey
Kuznetsov implemented in Fortran and provided on
[GitHub](https://github.com/Alexey-Kuznetsov-math/Riemann_zeta-and-Gamma).
They are usable in C++ as extern "C" functions.

Function | Description
--- | ---
czeta, czetaq | [Riemann zeta function ζ](https://en.wikipedia.org/wiki/Riemann_zeta_function) (double and quad precision implementations)
czetap, czetapq | Riemann zeta function ζ and its first derivative ζ′ (double and quad precision implementations)
clgammaq | Natural logarithm of the [gamma function Γ](https://en.wikipedia.org/wiki/Gamma_function) (quad precision only)
cpsiq | [Psi, alias digamma, function ψ](https://en.wikipedia.org/wiki/Digamma_function) = Γ′/Γ (quad precision only)

The names bound to the functions were chosen for consistency with the [C
standard library](https://en.cppreference.com/c/numeric) as far as that went
and then with the [GCC quad-precision math library
(libquadmath)](https://gcc.gnu.org/onlinedocs/gcc-16.1.0/libquadmath/)
thereafter.  See [notes below](#namconv) for an explanation of the naming
conventions.

The new functions are provided by libakzeta and declared in akzeta.h:

    typedef struct {
      double _Complex zeta;
      double _Complex zetap;
    } czetap_result;

    typedef struct {
      _Float128 _Complex zeta;
      _Float128 _Complex zetap;
    } czetapq_result;

       double _Complex czeta    (_Float128 _Complex s);
    _Float128 _Complex czetaq   (_Float128 _Complex s);
         czetap_result czetap   (_Float128 _Complex s);
        czetapq_result czetapq  (_Float128 _Complex s);
    _Float128 _Complex clgammaq (_Float128 _Complex z);
    _Float128 _Complex cpsiq    (_Float128 _Complex z);

Functions czeta and czetap use _Float128 internally, but they do shorter
calculations to produce fewer significant digits than czetaq and czetapq.
There would be no advantage in limiting their inputs to double precision.

A GNU environment (GNU make, gcc, gfortran, GNU libc, etc.) is assumed.  The
source code of included programs uses GNU extensions.

My additions are © 2026 David Flater but are provided under the same terms
as the originals (BSD 3-Clause License).

## Usable range and accuracy

|Im(s)| must not exceed 2.8e+19.  If it does, the library will just exit.

For |Im(s)| < 100, the expected accuracy of czetaq and czetapq is about 31
decimal digits.  As |Im(s)| increases, the number of good digits decreases to
approximately 32−log₁₀(|Im(s)|).

clgammaq is expected to maintain a relative accuracy of 10⁻³³.  cpsiq is
expected to maintain a relative accuracy of 10⁻³² except when z is close to
the negative real axis, where accuracy degrades.

Practically, accuracy stated in such terms does not apply in the vicinity of
the trivial zeros, the poles of ζ and Γ, or far into the negative half-plane
where 128-bit floats overflow (Re(s) < −2312.5).

The accuracy of czeta and czetap is capped at 15 digits since that is the
limit of the double precision data type.

## Algorithms

The zeta functions implement a new algorithm that is described in
[[preprint](https://arxiv.org/abs/2503.09519), [paywalled
DOI](https://doi.org/10.1016/j.cam.2026.117791)].  It combines the main sum
of
[Riemann–Siegel](https://en.wikipedia.org/wiki/Riemann%E2%80%93Siegel_formula)
with a simple approximation of the remainder term.

The zeta functions divide the domain as follows:

- For points close to the origin (abs(s) < 0.01), use a Taylor series.
- If Re(s) < ½, apply the reflection formula and continue with 1 − s.
- For |Im(s)| < 400, use whichever of
  [Euler–Maclaurin](https://en.wikipedia.org/wiki/Riemann_zeta_function#Numerical_algorithms)
  or a Dirichlet series has lower computational complexity for s.
- Otherwise, use whichever of Kuznetsov's algorithm or a Dirichlet series
  has lower computational complexity for s.

The gamma and psi functions implement a new alternative to [Lanczos
approximation](https://en.wikipedia.org/wiki/Lanczos_approximation) that is
described in [[preprint](https://arxiv.org/abs/2109.12061), [paywalled
DOI](https://doi.org/10.1016/j.cam.2022.114270)] and
[[preprint](https://arxiv.org/pdf/2508.19095),
[paywalled DOI](https://doi.org/10.1016/j.cam.2026.117756)].

[[Gourdon and Sebah
2003](http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.pdf)]
provides an introduction to previous algorithms for approximating zeta in
regions where a Dirichlet series summation does not suffice (Euler–Maclaurin,
Borwein, Riemann–Siegel, Odlyzko–Schönhage).  See also [[Borwein
1995](https://www.cecm.sfu.ca/~pborwein/PAPERS/P155.pdf)] for Borwein
Algorithms 2 and 3, which Gourdon and Sebah referred to as Proposition 1 and
2 under the convergence of alternating series method.

## Demo programs

- czetap:  Test czetap from the command line.
Usage: czetap sreal simag
- czetapq:  Test czetapq from the command line.
Usage: czetapq sreal simag
- clgammaq:  Test clgammaq from the command line.
Usage: clgammaq zreal zimag
- cpsiq:  Test cpsiq from the command line.
Usage: cpsiq zreal zimag
- test:  Run noninteractive test suite.
Usage: test

To build all demo programs normally, just type `make`.

To build all demo programs in C++ mode, type `make cxx`.  The C++ binaries
will be called czetap_cxx, czetapq_cxx, etc.

## Sample output

- [test_out.txt](test_out.txt):  Sample output of test.c.
- [test_out_upstream_Fortran.txt](test_out_upstream_Fortran.txt):
Sample output of the upstream test.f90 built according to the author's
instructions.
- [test_out_upstream_Python_Fortran.txt](test_out_upstream_Python_Fortran.txt):
Sample output of the upstream test.py built according to the author's
instructions.

## Notes

### <a name="namconv"> Naming conventions

Explaining how the names of the functions are consistent with the [C standard
library](https://en.cppreference.com/c/numeric) and the [GCC quad-precision
math library
(libquadmath)](https://gcc.gnu.org/onlinedocs/gcc-16.1.0/libquadmath/).

Point 0:  function names are ASCII; no Greek letters.

Point 1:  the real [log gamma
function](https://en.cppreference.com/c/numeric/math/lgamma) is called
lgamma:

          float lgammaf (float)
         double lgamma  (double)
    long double lgammal (long double)
      _Float128 lgammaq (_Float128)

The C standard library does not yet provide real or complex versions of
zeta, zeta prime, or psi.  Some names from other libraries:

Function | [GNU Scientific Library](https://www.gnu.org/software/gsl/doc/html/specfunc.html) | [C++ standard library](https://en.cppreference.com/cpp/numeric/special_functions) | [Boost Math](https://www.boost.org/doc/libs/latest/libs/math/doc/html/special.html) | [Python mpmath](https://mpmath.org)
--- | --- | --- | --- | ---
ζ     | gsl_sf_zeta    | riemann_zeta | zeta    | zeta
ζ′    |                |              |         | zeta (with derivative=1)
log Γ | gsl_sf_lngamma | lgamma       | lgamma  | loggamma
ψ     | gsl_sf_psi     |              | digamma | psi (with m=0) or digamma

Point 2:  the precision is indicated by an optional 1-letter suffix:

Suffix | C data type | Meaning on x86
--- | --- | ---
f | float | 32 bits (single precision)
  | double | 64 bits (double precision)
l | long double | 80 bits (extended precision)
q | _Float128 | 128 bits (quad precision)

One also finds suffixes of the form f<i>N</i> where <i>N</i> is the number of
bits in the floating-point type.  This competing naming convention comes from
[ISO/IEC TS 18661-3](https://open-std.org/JTC1/SC22/WG14/www/docs/n2601.pdf)
(Floating-point extensions for C, Part 3: Interchange and extended types) and
from [C++](https://en.cppreference.com/cpp/language/floating_literal).

Point 3:  the names of functions taking complex numbers are formed by
prefixing a c to the names of the corresponding real functions; e.g., [clogf,
clog, and clogl](https://en.cppreference.com/c/numeric/complex/clog) as the
complex analogs of [logf, log, and
logl](https://en.cppreference.com/c/numeric/math/log).

### Alternatives

Other free implementations of zeta for complex s:

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

Kuznetsov's ln_gamma is more complex than a "standard" Lanczos approximation,
for which there are many double-precision implementations using [various
widely-copied sets of
coefficients](https://www.mrob.com/pub/ries/lanczos-gamma.html).  But
apparently, the standard approximation is hard to implement at quad
precision:

- [Boost](https://www.boost.org/doc/libs/latest/libs/math/doc/html/math_toolkit/lanczos.html)
had to convert the Lanczos sum into rational form to obtain satisfactory
results for 128-bit floats.  "That means that the sum in rational form can be
evaluated without cancellation error, albeit with double the number of
coefficients for a given N."
- [libquadmath/math/lgammaq.c](https://github.com/gcc-mirror/gcc/blob/master/libquadmath/math/lgammaq.c)
partitions the domain into many segments, resulting in a much longer
function.  I have not found any notes or documentation for this
implementation.

### Floating point type traits

[Standard C library header float.h](https://en.cppreference.com/c/header/float)

lib/gcc/x86_64-pc-linux-gnu/17.0.0/include/float.h

DIG:  Number of decimal digits, q, such that any floating-point number with q
decimal digits can be rounded into a floating-point number with p radix b
digits and back again without change to the q decimal digits.  Number of
decimal digits that are guaranteed to be preserved in text →
float/double/long double → text roundtrip without change due to rounding or
overflow.

DECIMAL_DIG:  Number of decimal digits, n, such that any floating-point
number in the widest supported floating type with pmax radix b digits can be
rounded to a floating-point number with n decimal digits and back again
without change to the value.  Conversion from float/double/long double to
decimal with at least FLT_DECIMAL_DIG/DBL_DECIMAL_DIG/LDBL_DECIMAL_DIG digits
and back is the identity conversion:  this is the decimal precision required
to serialize/deserialize a floating-point value.

type | sizeof | prefix | DIG | DECIMAL_DIG
---: | ---: | ---: | ---: | ---:
   _Float16 |  2 |  FLT16 |  3 |  5
   _Float32 |  4 |  FLT32 |  6 |  9
      float |  4 |    FLT |  6 |  9
   _Float64 |  8 |  FLT64 | 15 | 17
     double |  8 |    DBL | 15 | 17
long double | 16 |   LDBL | 18 | 21
  _Float128 | 16 | FLT128 | 33 | 36
