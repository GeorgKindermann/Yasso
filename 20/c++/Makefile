#OP = -Ofast -march=native
OP = -O3 -march=native

WR = -Wall -Wextra

all:	y20c_subroutine.o simpleExample

y20c_subroutine.o:	y20c_subroutine.h y20c_subroutine.cc
	g++ -c $(OP) $(WR) y20c_subroutine.cc

simpleExample:	simpleExample.cc y20c_subroutine.o
	g++ $(OP) $(WR) -std=c++17 simpleExample.cc -osimpleExample y20c_subroutine.o
	strip simpleExample

clean:
	rm simpleExample y20c_subroutine.o
