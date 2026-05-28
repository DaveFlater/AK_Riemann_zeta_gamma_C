// © 2026 David Flater BSD-3-Clause

#include "common.h"

int main (int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: czetap sreal simag\n");
    exit(EX_USAGE);
  }
  init_stdout();
  const uint8_t namewidth = 5U;
  _Float128 sr, si;
  scan_f128(argv[1], &sr);
  scan_f128(argv[2], &si);
  _Float128 _Complex s = CMPLXF128(sr, si);
  print_f128_complex(s, L"s", namewidth);
  czetap_result ret = czetap(s);
  print_f64_complex(ret.zeta, L"ζ(s)", namewidth);
  print_f64_complex(ret.zetap, L"ζ′(s)", namewidth);
  return 0;
}
