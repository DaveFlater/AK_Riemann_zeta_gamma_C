// Functions provided by Riemann_zeta.f90

#pragma once

#ifdef __cplusplus
  extern "C" _Float128 _Complex Riemann_zeta (_Float128 _Complex s);
  extern "C" _Float128 _Complex ln_gamma     (_Float128 _Complex z);
#else
  extern _Float128 _Complex Riemann_zeta (_Float128 _Complex s);
  extern _Float128 _Complex ln_gamma     (_Float128 _Complex z);
#endif
