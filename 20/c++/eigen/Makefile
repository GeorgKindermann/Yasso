CXX = g++
#CXX = clang++

OP = -Ofast -march=native
#OP = -O3 -march=native
#OP = -O3 -ffast-math -march=native -fvectorize -fslp-vectorize

WR = -Wall -Wextra

STD = -std=c++17

all:	y20c_subroutine.o simpleExample

y20c_subroutine.o: y20c_subroutine.h y20c_subroutine.cc
	$(CXX) -c $(OP) $(WR) $(STD) -I /usr/include/eigen3 y20c_subroutine.cc

simpleExample: simpleExample.cc y20c_subroutine.o
	$(CXX) $(OP) $(WR) $(STD) simpleExample.cc -o simpleExample y20c_subroutine.o
	strip simpleExample

clean:
	rm -f simpleExample y20c_subroutine.o
