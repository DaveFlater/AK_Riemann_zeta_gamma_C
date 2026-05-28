.EXTRA_PREREQS = Makefile

CC = gcc
CXX = g++
FC = gfortran
AR = gcc-ar
OPTFLAGS = -O2 -flto
CFLAGS = -std=gnu23 $(OPTFLAGS) -Wall -Wextra
CXXFLAGS = -std=gnu++26 $(OPTFLAGS) -Wall -Wextra
FFLAGS = $(OPTFLAGS) -Wall -Wextra -Wno-tabs
ARFLAGS =
LDFLAGS = $(OPTFLAGS) -L. -s

NORMALBINS = clgammaq cpsiq czetap czetapq test
CXXBINS = clgammaq_cxx cpsiq_cxx czetap_cxx czetapq_cxx test_cxx

normal: $(NORMALBINS)
cxx: $(CXXBINS)
all: normal cxx

# Built sources
f64.f90: upstream_src/Riemann_zeta_module.f90 f64_sed.txt
	sed -f f64_sed.txt upstream_src/Riemann_zeta_module.f90 > $@
f128.f90: upstream_src/gamma_zeta_module.f90 f128_sed.txt
	sed -f f128_sed.txt upstream_src/gamma_zeta_module.f90 > $@

libakzeta.a: f64.o f128.o
	$(AR) rcsv $(ARFLAGS) $@ $^

clgammaq.o:       clgammaq.c   common.h akzeta.h
clgammaq_cxx.o:   clgammaq.c   common.h akzeta.h
cpsiq.o:          cpsiq.c   common.h akzeta.h
cpsiq_cxx.o:      cpsiq.c   common.h akzeta.h
czetap.o:         czetap.c     common.h akzeta.h
czetap_cxx.o:     czetap.c     common.h akzeta.h
czetapq.o:        czetapq.c    common.h akzeta.h
czetapq_cxx.o:    czetapq.c    common.h akzeta.h
test.o:           test.c common.h akzeta.h
test_cxx.o:       test.c common.h akzeta.h
# plot.o:           plot.c common.h akzeta.h

clgammaq:      clgammaq.o common.o libakzeta.a
	$(FC) $(LDFLAGS) -o $@ $< common.o -lakzeta
clgammaq_cxx:  clgammaq_cxx.o common_cxx.o libakzeta.a
	$(FC) $(LDFLAGS) -o $@ $< common_cxx.o -lakzeta
cpsiq:         cpsiq.o common.o libakzeta.a
	$(FC) $(LDFLAGS) -o $@ $< common.o -lakzeta
cpsiq_cxx:     cpsiq_cxx.o common_cxx.o libakzeta.a
	$(FC) $(LDFLAGS) -o $@ $< common_cxx.o -lakzeta
czetap:        czetap.o common.o libakzeta.a
	$(FC) $(LDFLAGS) -o $@ $< common.o -lakzeta
czetap_cxx:    czetap_cxx.o common_cxx.o libakzeta.a
	$(FC) $(LDFLAGS) -o $@ $< common_cxx.o -lakzeta
czetapq:       czetapq.o common.o libakzeta.a
	$(FC) $(LDFLAGS) -o $@ $< common.o -lakzeta
czetapq_cxx:   czetapq_cxx.o common_cxx.o libakzeta.a
	$(FC) $(LDFLAGS) -o $@ $< common_cxx.o -lakzeta
test:          test.o common.o libakzeta.a
	$(FC) $(LDFLAGS) -o $@ $< common.o -lakzeta
test_cxx:      test_cxx.o common_cxx.o libakzeta.a
	$(FC) $(LDFLAGS) -o $@ $< common_cxx.o -lakzeta
# plot:          plot.o common.o libakzeta.a
# 	$(FC) $(LDFLAGS) -o $@ $< common.o -lakzeta -lpng

%.o : %.c
	$(CC) $(CFLAGS) -c $<

%_cxx.o : %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

README.html: README.md
	python3 -m markdown -x tables README.md > README.html

clean:
	rm -f $(NORMALBINS) $(CXXBINS) *.o *.mod libakzeta.a f64.f90 f128.f90 README.html
