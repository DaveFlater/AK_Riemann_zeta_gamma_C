// test.f90 by Alexey Kuznetsov
// Manually converted to C by David Flater 2026-03-09
// © Alexey Kuznetsov and David Flater BSD-3-Clause

#include "common.h"
#include <time.h>

void test_zeros (wchar_t const desc[], _Float128 const zeros[],
                 uint8_t n_zeros) {
  wprintf(L"\n");
  wprintf(desc);
  wprintf(L"\n");
  _Float128 maxerr = 0;
  wchar_t const * const fmt = (zeros[0] < 1e6 ?
                                (zeros[0] < 100 ? L"ζ(½+%8.5lfi)"
                                                : L"ζ(½+%8.4lfi)")
                                : L"ζ(½+…%7.5lfi)");
  for (int8_t k=0; k<n_zeros; ++k) {
    wchar_t name[15];
    swprintf(name, 15, fmt, fmod((double)zeros[k],1e6));
    const _Float128 _Complex computed_val =
      Riemann_zeta(CMPLXF128(0.5f128, zeros[k]));
    const _Float128 err = cabsf128(computed_val);
    print_f128_complex(computed_val, name);
    if (err > maxerr) maxerr = err;
  }
  print_f128(maxerr, L"Max error");
}

int main ([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
  init_stdout();

  {
    wprintf(L"Test 1:  compute Γ(k) for integer values of k\n");
    const uint32_t g[10] = {1U, 1U, 2U, 6U, 24U, 120U, 720U, 5040U, 40320U,
                            362880U};
    _Float128 maxerr = 0;
    for (uint8_t k=1; k<=10; ++k) {
      wchar_t name[10];
      swprintf(name, 10, L"Γ(%u)", k);
      const _Float128 _Complex
        computed_val = cexpf128(ln_gamma(CMPLXF128(k, 0))),
        expected_val = CMPLXF128(g[k-1], 0);
      const _Float128 err = cabsf128(computed_val / expected_val - 1);
      print_f128_complex(computed_val, name);
      if (err > maxerr) maxerr = err;
    }
    print_f128(maxerr, L"Max rel error");
  }{
    wprintf(L"\nTest 2:  compute Γ(k-½) for integer values of k\n");
    const _Float128 sqrtpi = sqrtf128(M_PIf128),
                      g[5] = {sqrtpi, 0.5f128*sqrtpi, 0.75f128*sqrtpi,
                              1.875f128*sqrtpi, 6.5625f128*sqrtpi};
    _Float128 maxerr = 0;
    for (uint8_t k=1; k<=5; ++k) {
      wchar_t name[10];
      swprintf(name, 10, L"Γ(%u-½)", k);
      const _Float128 _Complex
        computed_val = cexpf128(ln_gamma(CMPLXF128(k-0.5, 0))),
        expected_val = CMPLXF128(g[k-1], 0);
      const _Float128 err = cabsf128(computed_val / expected_val - 1);
      print_f128_complex(computed_val, name);
      if (err > maxerr) maxerr = err;
    }
    print_f128(maxerr, L"Max rel error");
  }{
    wprintf(L"\nTest 3:  compute |log(Γ(z+1))−log(Γ(z))−log(z)| for 10000 random values of z\n"
            "randomly distributed in the square |Im(z)|<20, |Re(z)|<20\n");
    _Float128 maxerr = 0;
    _Float128 _Complex maxerr_z = CMPLXF128(0,0);
    srand48((long)time(NULL));
    for (uint16_t k=0; k<10000U; ++k) {
      const double U1 = drand48(), U2 = drand48(); // [0.0, 1.0)
      const _Float128 _Complex z = CMPLXF128(40*U1-20, 40*U2-20),
        computed_val = ln_gamma(z+1) - ln_gamma(z) - clogf128(z);
      const _Float128 err = cabsf128(computed_val);
      if (err > maxerr) {
        maxerr = err;
        maxerr_z = z;
      };
    }
    print_f128(maxerr, L"Max error");
    print_f128_complex(maxerr_z, L"Max error z");
  }{
    wprintf(L"\nTest 4:  compute ζ(s) for s=2,4,6,8,10\n");
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
      1.64493406684822643647241516664602507f128,
      1.08232323371113819151600369654116782f128,
      1.01734306198444913971451792979092046f128,
      1.00407735619794433937868523850865221f128,
      1.00099457512781808533714595890031882f128};
    _Float128 maxerr = 0;
    for (uint8_t k=1; k<=5; ++k) {
      wchar_t name[10];
      swprintf(name, 10, L"ζ(%u)", 2*k);
      const _Float128 _Complex
        computed_val = Riemann_zeta(CMPLXF128(2*k, 0)),
        expected_val = CMPLXF128(g[k-1], 0);
      const _Float128 err = cabsf128(computed_val - expected_val);
      print_f128_complex(computed_val, name);
      if (err > maxerr) maxerr = err;
    }
    print_f128(maxerr, L"Max error");
  }{
    wprintf(L"\nTest 5:  compute ζ(s) for s=0,-1,-3,-5,-7,-9\n");
    const _Float128 g[6] = {-0.5f128, -1.0f128/12, 1.0f128/120,
			    -1.0f128/252, 1.0f128/240, -1.0f128/132};
    _Float128 maxerr = 0;
    for (int8_t k=0; k<=5; ++k) {
      const int8_t s = (k ? 1-2*k : 0);
      wchar_t name[10];
      swprintf(name, 10, L"ζ(%d)", s);
      const _Float128 _Complex
        computed_val = Riemann_zeta(CMPLXF128(s, 0)),
        expected_val = CMPLXF128(g[k], 0);
      const _Float128 err = cabsf128(computed_val - expected_val);
      print_f128_complex(computed_val, name);
      if (err > maxerr) maxerr = err;
    }
    print_f128(maxerr, L"Max error");
  }

  // The zeros of ζ(s) used in the next three tests were downloaded from
  // https://www-users.cse.umn.edu/~odlyzko/zeta_tables/index.html
  {
    _Float128 const zeros[5] = {
      14.13472514173469379045725198356247027f128,
      21.02203963877155499262847959389690277f128,
      25.01085758014568876321379099256282181f128,
      30.42487612585951321031189753058409132f128,
      32.93506158773918969066236896407490348f128
    };
    test_zeros(L"Test 6:  compute ζ(s) for the first five nontrivial zeros",
               zeros, 5);
  }{
    _Float128 const zeros[5] = {
      201.2647519437037887330161334275481732f128,
      202.4935945141405342776866606378643158f128,
      204.1896718031045543307164383863136851f128,
      205.3946972021632860252123793906930909f128,
      207.9062588878062098615019679077536442f128
    };
    test_zeros(L"Test 7:  compute ζ(s) for the first five nontrivial zeros above t=200",
               zeros, 5);
  }{
    _Float128 const zeros[5] = {
      1000000.584097696f128,
      1000000.828343490f128,
      1000001.435265267f128,
      1000001.905648406f128,
      1000002.877617740f128
    };
    test_zeros(L"Test 8:  compute ζ(s) for the first five nontrivial zeros above t=10⁶\n"
      "Note that these zeros are correct only to 10⁻⁹, so the values should be of\n"
      "the same order of magnitude.",
               zeros, 5);
  }

  return 0;
}
