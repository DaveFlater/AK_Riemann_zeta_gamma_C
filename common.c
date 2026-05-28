// © 2026 David Flater BSD-3-Clause
// Common includes for tests and demos

/*
  Infrastructure for parsing and printing values.  It's rather cumbersome.

  1.  GNU libc does not have scanf or printf length modifiers that work with
  _Float128.  There used to be support for %Q if libquadmath was in use, but
  it went away:
  https://sourceware.org/pipermail/libc-alpha/2023-February/145106.html

  2.  For field widths to be calculated correctly when there are Greek
  letters in the string, you need to use wide chars.
*/

#include "common.h"

void init_stdout () {
  setlocale(LC_CTYPE, "en_US.UTF-8");
  if (fwide(stdout, 1) <= 0)
    fprintf(stderr, "fwide failed in init_stdout\n");
}

// strfromf128 is in stdlib.h but undocumented; see man strfromf.
// "%+.36g"  The + is unsupported and it aborts the program.
void format_f128 (_Float128 f, char buf[80]) {
  (void)strfromf128(buf, 79, "%.36g", f);
  buf[79] = '\0';
}

void print_f128 (_Float128 f, wchar_t const name[], uint8_t namewidth) {
  char fs[80];
  format_f128(f, fs);
  if (wprintf(L"%*ls = %s\n", namewidth, name, fs) < 1)
    fprintf(stderr, "wprintf failed in print_f128\n");
}

void print_f128_complex (_Float128 _Complex c, wchar_t const name[],
			 uint8_t namewidth) {
  char cr[80], ci[80];
  format_f128(Re(c), cr);
  format_f128(Im(c), ci);
  if (wprintf(L"%*ls = %s%s%si\n", namewidth, name, cr,
	      (Im(c) >= 0 ? "+" : ""), ci) < 1)
    fprintf(stderr, "wprintf failed in print_f128_complex\n");
}

void print_f64_complex (double _Complex c, wchar_t const name[],
			uint8_t namewidth) {
  if (wprintf(L"%*ls = %.17g%+.17gi\n", namewidth, name, Re(c), Im(c)) < 1)
    fprintf(stderr, "wprintf failed in print_f64_complex\n");
}

// strtof128 is in stdlib.h but undocumented; see man strtod.
void scan_f128 (char const buf[], _Float128 *f) {
  char *endptr;
  errno = 0;
  *f = strtof128(buf, &endptr);
  if (endptr == buf)
    fprintf(stderr, "Failed to convert to _Float128: %s\n", buf);
  if (errno == ERANGE)
    fprintf(stderr, "Range error (%sflow) on conversion of: %s\n",
	    (*f == HUGE_VAL_F128 || *f == -HUGE_VAL_F128 ? "over" : "under"),
	    buf);
  if (endptr == buf || errno == ERANGE) {
    print_f128(*f, L"Value", 5U);
    exit(EX_DATAERR);
  }
}
