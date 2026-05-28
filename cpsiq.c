// © 2026 David Flater BSD-3-Clause

#include "common.h"

int main (int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: cpsiq zreal zimag\n");
    exit(EX_USAGE);
  }
  init_stdout();
  _Float128 zr, zi;
  scan_f128(argv[1], &zr);
  scan_f128(argv[2], &zi);
  _Float128 _Complex s = CMPLXF128(zr, zi);
  print_f128_complex(s, L"z", 4U);
  _Float128 _Complex ret = cpsiq(s);
  print_f128_complex(ret, L"ψ(z)", 4U);
  return 0;
}
