// © 2026 David Flater BSD-3-Clause
// Common includes for zeta.c, gamma.c, and test.c.

/*
  Infrastructure for parsing and printing values.  It's rather cumbersome.

  1.  GNU libc does not have scanf or printf length modifiers that work with
  _Float128.  There used to be support for %Q if libquadmath was in use, but
  it went away:
  https://sourceware.org/pipermail/libc-alpha/2023-February/145106.html

  2.  For field widths to be calculated correctly when there are Greek
  letters in the string, you need to use wide chars.
*/

#pragma once

#include "Riemann_zeta.h"
#include <sysexits.h>

#ifdef __cplusplus

  #include <complex>
  #include <cerrno>
  #include <cinttypes>
  #include <clocale>
  #include <cmath>
  #include <cstdio>
  #include <cstdlib>
  #include <cwchar>
  extern "C" void init_stdout ();
  extern "C" void format_f128 (_Float128 f, char buf[80]);
  extern "C" void print_f128 (_Float128 f, wchar_t const name[]);
  extern "C" void print_f128_complex (_Float128 _Complex c,
				      wchar_t const name[]);
  extern "C" void scan_f128 (char const buf[], _Float128 *f);

  // This macro and declarations disappear in C++ mode.
  #define CMPLXF128(x, y) {(_Float128)(x), (_Float128)(y)}
  extern "C" _Float128 _Complex cexpf128 (_Float128 _Complex z);
  extern "C" _Float128 _Complex clogf128 (_Float128 _Complex z);
  extern "C" _Float128          cabsf128 (_Float128 _Complex z);

#else

  #define _GNU_SOURCE
  #include <complex.h>
  #include <errno.h>
  #include <inttypes.h>
  #include <locale.h>
  #include <math.h>
  #include <stdio.h>
  #include <stdlib.h>
  #include <wchar.h>
  extern void init_stdout ();
  extern void format_f128 (_Float128 f, char buf[80]);
  extern void print_f128 (_Float128 f, wchar_t const name[]);
  extern void print_f128_complex (_Float128 _Complex c,
				  wchar_t const name[]);
  extern void scan_f128 (char const buf[], _Float128 *f);

#endif
