#OP = -Ofast -march=native
OP = -O3 -march=native

WR = -Wall -Wextra

all:	yassofortran.o y20c_subroutine.o simpleExample

yassofortran.o:	yassofortran.f90
	gfortran -c $(OP) yassofortran.f90

y20c_subroutine.o:	y20c_subroutine.h y20c_subroutine.cc
	g++ -c $(OP) $(WR) y20c_subroutine.cc

simpleExample:	simpleExample.cc yassofortran.o y20c_subroutine.o
	g++ $(OP) $(WR) -std=c++17 simpleExample.cc -osimpleExample yassofortran.o y20c_subroutine.o -lgfortran -L/usr/local/lib64/ -Wl,-rpath -Wl,/usr/local/lib64/
	strip simpleExample

clean:
	rm yassofortran.o simpleExample yasso20.mod y20c_subroutine.o
