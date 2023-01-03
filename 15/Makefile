#OP = -Ofast -march=native
OP = -O3 -march=native

WR = -Wall -Wextra

all:	y15_subroutine.o y15c_subroutine.o simpleExample simplePerformance

<<<<<<< HEAD:15/c++/Makefile
y15_subroutine.o:	../fortran/y15_subroutine.f90
	gfortran -c $(OP) ../fortran/y15_subroutine.f90
=======
y15_subroutine.o:	y15_subroutine.f90
	gfortran -c $(OP) y15_subroutine.f90
>>>>>>> 140064075045444e2d0fe1d14ac8a9b3060cbc11:15/Makefile

y15c_subroutine.o:	y15c_subroutine.h y15c_subroutine.cc
	g++ -c $(OP) $(WR) -std=c++17 y15c_subroutine.cc

simpleExample:	simpleExample.cc y15_subroutine.o y15c_subroutine.o
	g++ $(OP) $(WR) -std=c++17 simpleExample.cc -osimpleExample y15_subroutine.o y15c_subroutine.o -lgfortran -L/usr/local/lib64/ -Wl,-rpath -Wl,/usr/local/lib64/
	strip simpleExample

simplePerformance: simplePerformance.cc y15_subroutine.o y15c_subroutine.o
	g++ $(OP) $(WR) -std=c++17 simplePerformance.cc -osimplePerformance y15_subroutine.o y15c_subroutine.o -lgfortran -L/usr/local/lib64/ -Wl,-rpath -Wl,/usr/local/lib64/ -Wno-unused-variable
	strip simplePerformance

clean:
	rm y15_subroutine.o y15c_subroutine.o simpleExample simplePerformance yasso.mod
