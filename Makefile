.EXTRA_PREREQS = Makefile

# Normal, build as C
CC = gcc
CFLAGS = -std=gnu23 -O2 -Wall -Wextra
normal:  zeta gamma test
zeta.o:  zeta.c  Riemann_zeta.h common.h
gamma.o: gamma.c Riemann_zeta.h common.h
test.o:  test.c  Riemann_zeta.h common.h
zeta:    zeta.o  Riemann_zeta.o common.o
gamma:   gamma.o Riemann_zeta.o common.o
test:    test.o  Riemann_zeta.o common.o

# Optional, test build as C++
CXX = g++
CXXFLAGS = -std=gnu++26 -O2 -Wall -Wextra
cxx:         zeta_cxx gamma_cxx test_cxx
zeta_cxx.o:  zeta.c      Riemann_zeta.h common.h
gamma_cxx.o: gamma.c     Riemann_zeta.h common.h
test_cxx.o:  test.c      Riemann_zeta.h common.h
zeta_cxx:    zeta_cxx.o  Riemann_zeta.o common_cxx.o
gamma_cxx:   gamma_cxx.o Riemann_zeta.o common_cxx.o
test_cxx:    test_cxx.o  Riemann_zeta.o common_cxx.o

all: normal cxx

%.o : %.c
	$(CC) $(CFLAGS) -c $<

%_cxx.o : %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

% : %.o
	gfortran -o $@ $^

Riemann_zeta.o: Riemann_zeta.f90
	gfortran -O2 -Wall -Wextra -Wno-real-q-constant -c Riemann_zeta.f90

README.html: README.md
	python3 -m markdown -x tables README.md > README.html

clean:
	rm -f zeta gamma test zeta_cxx gamma_cxx test_cxx Riemann_zeta.o zeta.o gamma.o test.o common.o zeta_cxx.o gamma_cxx.o test_cxx.o common_cxx.o riemann_zeta_module.mod README.html

# Quick and dirty, does not put the files in a subdirectory
dist:
	tar cvJf Riemann_zeta.tar.xz Makefile *.f90 *.c *.h *.txt *.md
