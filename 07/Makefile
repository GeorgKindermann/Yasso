#OP = -Ofast -march=native
OP = -O3 -march=native

WR = -Wall -Wextra

all:	y07_subroutine.o simpleExample

y07_subroutine.o:	y07_subroutine.f90
	gfortran -c $(OP) y07_subroutine.f90

simpleExample:	simpleExample.cc y07_subroutine.o
	g++ $(OP) $(WR) -std=c++17 simpleExample.cc -osimpleExample y07_subroutine.o -lgfortran -L/usr/local/lib64/ -Wl,-rpath -Wl,/usr/local/lib64/
	strip simpleExample

clean:
	rm y07_subroutine.o simpleExample yasso07.mod
