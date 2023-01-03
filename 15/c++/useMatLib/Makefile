OP = -Ofast -march=native
#OP = -O3 -march=native

#WR = -Wall -Wextra

MKLROOT = /opt/intel/mkl/

all:	 y15cArma_subroutine.o simplePerformanceArma y15cEigen_subroutine.o simplePerformanceEigen y15cAtlas_subroutine.o simplePerformanceAtlas y15cCuda_subroutine.o simplePerformanceCuda y15cMkl_subroutine.o simplePerformanceMkl

y15cArma_subroutine.o:	y15cArma_subroutine.h  y15cArma_subroutine.cc
	g++ -c $(OP) $(WR) -std=c++17  y15cArma_subroutine.cc

y15cEigen_subroutine.o:	y15cEigen_subroutine.h  y15cEigen_subroutine.cc
	g++ -c $(OP) $(WR) -std=c++17  y15cEigen_subroutine.cc -I/usr/include/eigen3

y15cAtlas_subroutine.o:	y15cAtlas_subroutine.h  y15cAtlas_subroutine.cc
	g++ -c $(OP) $(WR) -std=c++17  y15cAtlas_subroutine.cc

y15cCuda_subroutine.o:	y15cCuda_subroutine.h  y15cCuda_subroutine.cc
	g++ -c $(OP) $(WR) -std=c++17  y15cCuda_subroutine.cc

y15cMkl_subroutine.o:	y15cMkl_subroutine.h  y15cMkl_subroutine.cc
	g++ -c $(OP) $(WR) -std=c++17 -DMKL_ILP64 -m64 -I${MKLROOT}/include y15cMkl_subroutine.cc

simplePerformanceArma: simplePerformanceArma.cc  y15cArma_subroutine.o
	g++ $(OP) $(WR) -std=c++17 simplePerformanceArma.cc -osimplePerformanceArma  y15cArma_subroutine.o -L/usr/local/lib64/ -Wl,-rpath -Wl,/usr/local/lib64/ -Wno-unused-variable -larmadillo
	strip simplePerformanceArma

simplePerformanceEigen: simplePerformanceEigen.cc  y15cEigen_subroutine.o
	g++ $(OP) $(WR) -std=c++17 simplePerformanceEigen.cc -osimplePerformanceEigen  y15cEigen_subroutine.o -L/usr/local/lib64/ -Wl,-rpath -Wl,/usr/local/lib64/ -Wno-unused-variable -I/usr/include/eigen3
	strip simplePerformanceEigen

simplePerformanceAtlas: simplePerformanceAtlas.cc  y15cAtlas_subroutine.o
	g++ $(OP) $(WR) -std=c++17 simplePerformanceAtlas.cc -osimplePerformanceAtlas  y15cAtlas_subroutine.o -L/usr/local/lib64/ -Wl,-rpath -Wl,/usr/local/lib64/ -Wno-unused-variable -L/usr/lib/atlas-base/atlas -lblas
	strip simplePerformanceAtlas

simplePerformanceCuda: simplePerformanceCuda.cc  y15cCuda_subroutine.o
	g++ $(OP) $(WR) -std=c++17 simplePerformanceCuda.cc -osimplePerformanceCuda  y15cCuda_subroutine.o -L/usr/local/lib64/ -Wl,-rpath -Wl,/usr/local/lib64/ -Wno-unused-variable -L/usr/lib/x86_64-linux-gnu/ -lcublas
	strip simplePerformanceCuda

simplePerformanceMkl: simplePerformanceMkl.cc  y15cMkl_subroutine.o
	g++ $(OP) $(WR) -std=c++17 simplePerformanceMkl.cc -osimplePerformanceMkl  y15cMkl_subroutine.o -L/usr/local/lib64/ -Wl,-rpath -Wl,/usr/local/lib64/ -Wno-unused-variable -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
	strip simplePerformanceMkl

clean:
	rm  y15cArma_subroutine.o simplePerformanceArma y15cEigen_subroutine.o simplePerformanceEigen y15cAtlas_subroutine.o simplePerformanceAtlas y15cCuda_subroutine.o simplePerformanceCuda  y15cMkl_subroutine.o simplePerformanceMkl
