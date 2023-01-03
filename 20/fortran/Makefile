#OP = -Ofast -march=native
OP = -O3 -march=native

WR = -Wall -Wextra

all:	yassofortran.o simpleExample

yassofortran.o:	yassofortran.f90
	gfortran -c $(OP) yassofortran.f90

simpleExample:	simpleExample.cc yassofortran.o
	g++ $(OP) $(WR) simpleExample.cc -osimpleExample yassofortran.o -lgfortran -L/usr/local/lib64/ -Wl,-rpath -Wl,/usr/local/lib64/
	strip simpleExample

clean:
	rm yassofortran.o simpleExample yasso20.mod
