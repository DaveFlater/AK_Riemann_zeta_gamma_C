// Fortran/test.f90 and Python_Fortran/test.py by Alexey Kuznetsov, manually
// converted and merged by David Flater
// © Alexey Kuznetsov and David Flater BSD-3-Clause

#include "common.h"
#include <time.h>
#include <assert.h>

void fmt_zero_name (_Float128 zero, uint16_t numwidth,
		    wchar_t name[], uint16_t namelen) {
  char numbuf[numwidth+1];
  // man strfromf says:  The terminating null byte ('\0') is written if and
  // only if n is sufficiently large, otherwise the written string is
  // truncated at n characters.  But strfromf128 is putting a null instead of
  // the fifth decimal when n == numwidth.
  (void)strfromf128(numbuf, numwidth+1, "%.5f", zero);
  numbuf[numwidth] = '\0';
  swprintf(name, namelen, L"ζ(½+%si)", numbuf);
}

void test_zeta_zeros (wchar_t const desc[],
		      _Float128 const zeros[],
		      uint8_t n_zeros) {
  wprintf(L"\n");
  wprintf(desc);
  wprintf(L"\n");
  _Float128 maxerr = 0;
  assert(zeros[0] > 10);
  const uint16_t ldigs = (uint16_t)log10f128(zeros[0])+1U,
              numwidth = ldigs + 6U,
             namewidth = numwidth + 6U,
               namelen = namewidth + 1U;
  wchar_t name[namelen];
  wprintf(L"czetaq:\n");
  for (int8_t k=0; k<n_zeros; ++k) {
    fmt_zero_name(zeros[k], numwidth, name, namelen);
    const _Float128 _Complex computed_val = czetaq(CMPLXF128(0.5q, zeros[k]));
    const _Float128 err = cabsf128(computed_val);
    print_f128_complex(computed_val, name, namewidth);
    if (err > maxerr) maxerr = err;
  }
  print_f128(maxerr, L"Max error", namewidth);
  maxerr = 0;
  wprintf(L"czeta:\n");
  for (int8_t k=0; k<n_zeros; ++k) {
    fmt_zero_name(zeros[k], numwidth, name, namelen);
    const double _Complex computed_val = czeta(CMPLXF128(0.5q, zeros[k]));
    const _Float128 err = cabsf128(computed_val);
    print_f64_complex(computed_val, name, namewidth);
    if (err > maxerr) maxerr = err;
  }
  print_f128(maxerr, L"Max error", namewidth);
}

void test_zeta_ints (wchar_t const desc[],
		     int const ints[],
		     _Float128 const expected[],
		     uint8_t n_ints) {
  wprintf(L"\n");
  wprintf(desc);
  wprintf(L"\n");
  _Float128 maxerr = 0;
  // The ints need 2 chars, -9 to 10
  wchar_t const fmt[] = L"ζ(%d)";
  wprintf(L"czetaq:\n");
  for (int8_t k=0; k<n_ints; ++k) {
    wchar_t name[6];
    swprintf(name, 6, fmt, ints[k]);
    const _Float128 _Complex
      computed_val = czetaq(CMPLXF128(ints[k], 0)),
      expected_val = CMPLXF128(expected[k], 0);
    const _Float128 err = cabsf128(computed_val - expected_val);
    print_f128_complex(computed_val, name, 9U);
    if (err > maxerr) maxerr = err;
  }
  print_f128(maxerr, L"Max error", 9U);
  maxerr = 0;
  wprintf(L"czeta:\n");
  for (int8_t k=0; k<n_ints; ++k) {
    wchar_t name[6];
    swprintf(name, 6, fmt, ints[k]);
    const double _Complex computed_val = czeta(CMPLXF128(ints[k], 0));
    const _Float128 _Complex expected_val = CMPLXF128(expected[k], 0);
    const _Float128 err = cabsf128(computed_val - expected_val);
    print_f64_complex(computed_val, name, 9U);
    if (err > maxerr) maxerr = err;
  }
  print_f128(maxerr, L"Max error", 9U);
}

void test_Newton (wchar_t const desc[],
		  _Float128 _Complex s,
		  _Float128 _Complex z) {
  wprintf(L"\n");
  wprintf(desc);
  wprintf(L"\n");
  {
    wprintf(L"czetapq:\n");
    _Float128 _Complex sq = s;
    print_f128_complex(sq, L"s", 5U);
    for (uint8_t k=0; k<5; ++k) {
      czetapq_result h = czetapq(sq);
      sq -= h.zeta/h.zetap;
      print_f128_complex(sq, L"s", 5U);
    }
    print_f128(cabsf128(sq-z), L"Error", 5U);
    print_f128_complex(czetaq(sq), L"ζ(s)", 5U);
  }{
    wprintf(L"czetap:\n");
    double _Complex sd = s;
    print_f64_complex(sd, L"s", 5U);
    for (uint8_t k=0; k<4; ++k) {
      czetap_result h = czetap(sd);
      sd -= h.zeta/h.zetap;
      print_f64_complex(sd, L"s", 5U);
    }
    print_f128(cabsf128(sd-z), L"Error", 5U);
    print_f64_complex(czeta(sd), L"ζ(s)", 5U);
  }
}

int main ([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
  init_stdout();
  wprintf(L"\
This test suite is a merger of tests from Fortran/test.f90 and\n\
Python_Fortran/test.py upstream.  Test numbers x/y refer to the test numbers\n\
in Fortran/test.f90 and Python_Fortran/test.py respectively.  Significant\n\
changes:\n\
- Used more precise expected results for Tests 10/4 and -/6.\n\
- Used the same (quad precision) input for both double and quad precision\n\
  implementations.\n\
- Calculated error in quad precision for both double and quad precision\n\
  implementations.\n");

  // 5 tests pertaining to lgamma and psi, quad precision only

  {
    wprintf(L"\nTest 1/-:  compute Γ(k) for integer values of k\n");
    const uint32_t g[10] = {1U, 1U, 2U, 6U, 24U, 120U, 720U, 5040U, 40320U,
                            362880U};
    _Float128 maxerr = 0;
    for (uint8_t k=1; k<=10; ++k) {
      wchar_t name[10];
      swprintf(name, 10, L"Γ(%u)", k);
      const _Float128 _Complex
        computed_val = cexpf128(clgammaq(CMPLXF128(k, 0))),
        expected_val = CMPLXF128(g[k-1], 0);
      const _Float128 err = cabsf128(computed_val / expected_val - 1);
      print_f128_complex(computed_val, name, 13U);
      if (err > maxerr) maxerr = err;
    }
    print_f128(maxerr, L"Max rel error", 13U);
  }{
    wprintf(L"\nTest 2/-:  compute Γ(k-½) for integer values of k\n");
    const _Float128 sqrtpi = sqrtf128(M_PIf128),
                      g[5] = {sqrtpi, 0.5q*sqrtpi, 0.75q*sqrtpi,
                              1.875q*sqrtpi, 6.5625q*sqrtpi};
    _Float128 maxerr = 0;
    for (uint8_t k=1; k<=5; ++k) {
      wchar_t name[10];
      swprintf(name, 10, L"Γ(%u-½)", k);
      const _Float128 _Complex
        computed_val = cexpf128(clgammaq(CMPLXF128(k-0.5, 0))),
        expected_val = CMPLXF128(g[k-1], 0);
      const _Float128 err = cabsf128(computed_val / expected_val - 1);
      print_f128_complex(computed_val, name, 13U);
      if (err > maxerr) maxerr = err;
    }
    print_f128(maxerr, L"Max rel error", 13U);
  }{
    wprintf(L"\nTest 3/-:  compute |log(Γ(z+1))−log(Γ(z))−log(z)| for 10000 random values of z\n"
            "randomly distributed in the square |Im(z)|<20, |Re(z)|<20\n");
    _Float128 maxerr = 0;
    _Float128 _Complex maxerr_z = CMPLXF128(0,0);
    srand48((long)time(NULL));
    for (uint16_t k=0; k<10000U; ++k) {
      const double U1 = drand48(), U2 = drand48(); // [0.0, 1.0)
      const _Float128 _Complex z = CMPLXF128(40*U1-20, 40*U2-20),
        computed_val = clgammaq(z+1) - clgammaq(z) - clogf128(z);
      const _Float128 err = cabsf128(computed_val);
      if (err > maxerr) {
        maxerr = err;
        maxerr_z = z;
      };
    }
    print_f128(maxerr, L"Max error", 11U);
    print_f128_complex(maxerr_z, L"Max error z", 11U);
  }{
    wprintf(L"\nTest 4/-:  check that psi(1)=-0.57721566490153286060651209008240243 (Euler-Mascheroni constant)\n");
    const _Float128 _Complex z = CMPLXF128(1,0),
      ret = cpsiq(z),
      err = ret + 0.57721566490153286060651209008240243q;
    print_f128_complex(ret, L"ψ(1)", 5U);
    print_f128_complex(err, L"Error", 5U);
  }{
    wprintf(L"\nTest 5/-:  compute |ψ(z+1)−ψ(z)−1/z| for 10000 random values of z\n"
            "randomly distributed in the square |Im(z)|<20, |Re(z)|<20\n");
    _Float128 maxerr = 0;
    _Float128 _Complex maxerr_z = CMPLXF128(0,0);
    srand48((long)time(NULL));
    for (uint16_t k=0; k<10000U; ++k) {
      const double U1 = drand48(), U2 = drand48(); // [0.0, 1.0)
      const _Float128 _Complex z = CMPLXF128(40*U1-20, 40*U2-20),
        computed_val = cpsiq(z+1) - cpsiq(z) - 1.0q/z;
      const _Float128 err = cabsf128(computed_val);
      if (err > maxerr) {
        maxerr = err;
        maxerr_z = z;
      };
    }
    print_f128(maxerr, L"Max error", 11U);
    print_f128_complex(maxerr_z, L"Max error z", 11U);
  }

  // Zeta tests, double and quad precision

  {
    /*
      powf128(M_PIf128, 10)/93555 evaluates to a slightly different value
      when test.c is built as C++ rather than C:
        C    1.00099457512781808533714595890031882
        C++  1.00099457512781808533714595890031863
      This causes a difference in max error even though the results from
      Riemann_zeta are identical.  To avoid that distraction, put in literal
      values.
      const _Float128 g[5] = {powf128(M_PIf128, 2)/6,
        powf128(M_PIf128, 4)/90,   powf128(M_PIf128, 6)/945,
        powf128(M_PIf128, 8)/9450, powf128(M_PIf128, 10)/93555};
    */
    const _Float128 g[5] = {
      1.64493406684822643647241516664602507q,
      1.08232323371113819151600369654116782q,
      1.01734306198444913971451792979092046q,
      1.00407735619794433937868523850865221q,
      1.00099457512781808533714595890031882q};
    const int ints[5] = {2, 4, 6, 8, 10};
    test_zeta_ints(L"Test 6/1:  compute ζ(s) for s=2,4,6,8,10",
		   ints, g, 5);
  }{
    const _Float128 g[6] = {-0.5q, -1.0q/12, 1.0q/120,
			    -1.0q/252, 1.0q/240, -1.0q/132};
    const int ints[6] = {0, -1, -3, -5, -7, -9};
    test_zeta_ints(L"Test 7/2:  compute ζ(s) for s=0,-1,-3,-5,-7,-9",
		   ints, g, 6);
  }

  // The zeros of ζ(s) used in the next two tests were downloaded from
  // https://www-users.cse.umn.edu/~odlyzko/zeta_tables/index.html
  {
    _Float128 const zeros[5] = {
      14.13472514173469379045725198356247027q,
      21.02203963877155499262847959389690277q,
      25.01085758014568876321379099256282181q,
      30.42487612585951321031189753058409132q,
      32.93506158773918969066236896407490348q
    };
    test_zeta_zeros(L"Test 8/3:  compute ζ(s) for the first five nontrivial zeros",
		    zeros, 5);
  }{
    _Float128 const zeros[5] = {
      201.2647519437037887330161334275481732q,
      202.4935945141405342776866606378643158q,
      204.1896718031045543307164383863136851q,
      205.3946972021632860252123793906930909q,
      207.9062588878062098615019679077536442q
    };
    test_zeta_zeros(L"Test 9/-:  compute ζ(s) for the first five nontrivial zeros above t=200",
		    zeros, 5);
  }
  // The zeros of zeta(s) used in the next two tests were downloaded from
  // https://www.lmfdb.org/zeros/zeta/
  {
    _Float128 const zeros[5] = {
      1000000.5840976963450700185233799609033q,
      1000000.8283434909528317159888785155369q,
      1000001.4352652671578816544819345694547q,
      1000001.9056484055705747742072156881134q,
      1000002.8776177398369455820195356920489q
    };
    test_zeta_zeros(L"Test 10/4:  compute ζ(s) for the first five nontrivial zeros above t=10⁶",
		    zeros, 5);
  }{
    _Float128 const zeros[5] = {
      1000000000.1156508900208481613883400760674q,
      1000000000.4340268958947434089036006725812q,
      1000000000.5303428567293815535254541945858q,
      1000000001.0190201774926669971382884553198q,
      1000000001.2937394927044076976924485252179q
    };
    test_zeta_zeros(L"Test -/5:  compute ζ(s) for the first five nontrivial zeros above t=10⁹",
		    zeros, 5);
  }

  {
    wprintf(L"\nTest -/6:  check 10 values of ζ\n");
    _Float128 maxerr = 0;
    const _Float128 _Complex s[10] = {
      CMPLXF128(-15.5q, 1.0q),
      CMPLXF128(0.5q, 20.0q),
      CMPLXF128(1.5q, -150.0q),
      CMPLXF128(0.0q, 500.0q),
      CMPLXF128(10.0q, 3000.0q),
      CMPLXF128(1.0q, -100000.0q),
      CMPLXF128(-0.5q, 5000000.0q),
      CMPLXF128(0.5q, 20000000.0q),
      CMPLXF128(0.75q, 1.0e+10q),
      CMPLXF128(0.4q, -1.0e+12q)
    };
    // Values from mpmath
    const _Float128 _Complex expected[10] = {
      CMPLXF128(1.60720992801538893636307075657527285461q, -0.314611728846899873791626296388441136793q),
      CMPLXF128(0.429913860437843372157739670624503456822q, -1.06429144308058911272739519306893847477q),
      CMPLXF128(0.643164286756546270307813834415938562918q, 0.138652708369380808021260267666799356449q),
      CMPLXF128(-9.82537217473325310405557723018451727036q, -2.18751366531488463083395942644030462596q),
      CMPLXF128(1.00091953780491926683453811452170842094q, 0.000287463103243799701078558890776821564166q),
      CMPLXF128(1.61812212284693679656775905790087726316q, -1.07044104147062368662603474543967113456q),
      CMPLXF128(1537232.27642331701352595512745291372431q, -590208.990876205560892555120517487148543q),
      CMPLXF128(0.28168940385795717592746083279566232929q, 2.35068498897859760614983495048006459377q),
      CMPLXF128(0.315056654699321886152305556563372986923q, 0.56144654896737661471557069874254476698q),
      CMPLXF128(14.3190966513808078283482956404477655357q, 17.8583347812398939685242501160951562115q)
    };
    wprintf(L"czetaq:\n");
    for (uint8_t i=0; i<10; ++i) {
      wchar_t name[80];
      swprintf(name, 80, L"ζ(%.2f%+.0fi)", (double)Re(s[i]), (double)Im(s[i]));
      const _Float128 _Complex computed_val = czetaq(s[i]);
      const _Float128 err = cabsf128(computed_val - expected[i]);
      print_f128_complex(computed_val, name, 22U);
      print_f128(err, L"Error", 22U);
      if (err > maxerr) maxerr = err;
    }
    print_f128(maxerr, L"Max error", 22U);
    maxerr = 0;
    wprintf(L"czeta:\n");
    for (uint8_t i=0; i<10; ++i) {
      wchar_t name[80];
      swprintf(name, 80, L"ζ(%.2f%+.0fi)", (double)Re(s[i]), (double)Im(s[i]));
      const double _Complex computed_val = czeta(s[i]);
      const _Float128 err = cabsf128(computed_val - expected[i]);
      print_f64_complex(computed_val, name, 22U);
      print_f128(err, L"Error", 22U);
      if (err > maxerr) maxerr = err;
    }
    print_f128(maxerr, L"Max error", 22U);
  }

  test_Newton(L"Test 11/7:  find the first non-trivial zero with Newton's method",
	      CMPLXF128(0.5q, 14.1q),
	      CMPLXF128(0.5q, 14.13472514173469379045725198356247027q));
  // Zero number 203
  test_Newton(L"Test 12/8:  find the zero near ½+401.8i with Newton's method",
	      CMPLXF128(0.5q, 401.8q),
	      CMPLXF128(0.5q, 401.83922860053321653991130437828723072q));

  // Skip Python_Fortran test 9, which is a plot of the Riemann-Siegel
  // function Z(t)

  return 0;
}
