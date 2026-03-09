// © 2026 David Flater BSD-3-Clause

#include "common.h"

int main (int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: zeta sreal simag\n");
    exit(EX_USAGE);
  }
  init_stdout();
  _Float128 sr, si;
  scan_f128(argv[1], &sr);
  scan_f128(argv[2], &si);
  _Float128 _Complex s = CMPLXF128(sr, si);
  print_f128_complex(s, L"s");
  _Float128 _Complex ret = Riemann_zeta(s);
  print_f128_complex(ret, L"ζ(s)");
  return 0;
}
