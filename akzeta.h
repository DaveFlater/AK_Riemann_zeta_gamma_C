// Functions provided by libakzeta

#pragma once

#ifdef __cplusplus
  #define LINKAGE extern "C"
  struct czetap_result {
    double _Complex zeta;
    double _Complex zetap;
  };
  struct czetapq_result {
    _Float128 _Complex zeta;
    _Float128 _Complex zetap;
  };
#else
  #define LINKAGE extern
  typedef struct {
    double _Complex zeta;
    double _Complex zetap;
  } czetap_result;
  typedef struct {
    _Float128 _Complex zeta;
    _Float128 _Complex zetap;
  } czetapq_result;
#endif

LINKAGE    double _Complex czeta    (_Float128 _Complex s);
LINKAGE _Float128 _Complex czetaq   (_Float128 _Complex s);
LINKAGE      czetap_result czetap   (_Float128 _Complex s);
LINKAGE     czetapq_result czetapq  (_Float128 _Complex s);
LINKAGE _Float128 _Complex clgammaq (_Float128 _Complex z);
LINKAGE _Float128 _Complex cpsiq    (_Float128 _Complex z);

#undef LINKAGE
