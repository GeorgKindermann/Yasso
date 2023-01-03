#include <cmath>
#include <complex>

#include "y15c_subroutine.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

  //https://github.com/ngocson2vn/MPM/blob/master/547_546-parallel-GIMP/0main/fit.cpp
  //http://www.sci.utah.edu/~wallstedt/LU.htm
  //Philip Wallstedt
  //Decompose
  template<const int d, typename T>inline void Crout(const std::array<T,d*d>& a, std::array<T,d*d>& A){
    for(int k=0;k<d;++k){
      for(int i=k;i<d;++i){
	T sum=0.;
	for(int p=0;p<k;++p)sum+=A[i*d+p]*A[p*d+k];
	A[i*d+k]=a[i*d+k]-sum;
      }
      for(int j=k+1;j<d;++j){
	T sum=0.;
	for(int p=0;p<k;++p)sum+=A[k*d+p]*A[p*d+j];
	A[k*d+j]=(a[k*d+j]-sum)/A[k*d+k];
      }
    }
  }
  template<const int d, typename T>inline void solveCrout(const std::array<T,d*d>& A, const std::array<T,d>& b, std::array<T,d>& x){
    std::array<T,d> y;
    for(int i=0;i<d;++i){
      T sum=0.;
      for(int k=0;k<i;++k)sum+=A[i*d+k]*y[k];
      y[i]=(b[i]-sum)/A[i*d+i];
    }
    for(int i=d-1;i>=0;--i){
      T sum=0.;
      for(int k=i+1;k<d;++k)sum+=A[i*d+k]*x[k];
      x[i]=(y[i]-sum);
    }
  }

  //Matrixmultiplication
  template<int M, int N, int P, typename T> void MATMUL(const std::array<T,M*N> &A, const std::array<T,N*P> &B, std::array<T,M*P> &C) {
    std::array<T,M*P> TMP;
    std::array<T,M*P>& RES = (B.data() == C.data() || A.data() == C.data()) ? TMP : C;
    for (int p=0; p<P; ++p) {
      for (int m=0; m<M; ++m) {
	T sum = A[m*N] * B[p];
	for (int n=1; n<N; ++n) {
	  sum += A[m*N+n] * B[n*P+p];
	}
	RES[m*P+p] = sum;
      }
    }
    if(&RES != &C) {for (int i=0; i<M*P; ++i) {C[i] = RES[i];}}
  }

  //Functions for solving the diff. equation, adapted for the Yasso case
  template<int M, typename T> T matrixnorm(const std::array<T, M*M>& A) {
    //returns elementwise (i.e. Frobenius) norm of a square matrix
    T b = 0.;
    for(int i=0; i<M*M; ++i) {b += A[i]*A[i];}
    return(sqrt(b));
  }
  template<int M, typename T> void matrixexp(const std::array<T, M*M>& A, std::array<T, M*M>& B, const size_t& q) {
    //Approximated matrix exponential using Taylor series with scaling & squaring
    //Accurate enough for the Yasso case
    //int q = 11; // #terms in Taylor  gk: possibility to reduce to speed up
    //expect B is initialized:  for(int i=0; i<25; ++i) {B[i] = 0.;}
    for(int i=0; i<M; ++i) {B[i*(M+1)] = 1.;}
    T normiter = 2.; // Amount of scaling & squaring
    int j=2;
    T p = matrixnorm<M>(A);
    while(p>normiter) {
      normiter *= 2.;
      ++j;
    }
    std::array<T, M*M> C,D;
    for(int i=0; i<M*M; ++i) {C[i] = A[i]/normiter;} //scale
    for(int i=0; i<M*M; ++i) {B[i] += C[i];}
    for(int i=0; i<M*M; ++i) {D[i] = C[i];}
    for(size_t i=2; i<q; ++i) { //compute Taylor expansion
      MATMUL<M,M,M>(C,D,D);
      for(int j=0; j<M*M; ++j) {D[j] /= T(i);}
      for(int i=0; i<M*M; ++i) {B[i] += D[i];}
    }
    for(int i=1; i<j; ++i) { //square
      MATMUL<M,M,M>(B,B,B);
    }
  }

}

namespace yasso {

  yasso15::yasso15() : A{0.} {
    setTheta({4.8971473e-01,4.9138734e+00,2.4197346e-01,9.4876416e-02,4.3628932e-01,2.4997402e-01,9.1512685e-01,9.9258227e-01,8.3853738e-02,1.1476783e-02,6.0831497e-04,4.7612821e-04,6.6037729e-02,7.7134168e-04,1.0401742e-01,6.4880756e-01  -1.5487177e-01  -1.9568024e-02  -9.1717130e-01  -4.0359430e-04  -1.6707272e-04,9.0598047e-02  -2.1440956e-04,4.8772465e-02  -7.9136021e-05,3.5185492e-02  -2.0899057e-04  -1.8089202e+00  -1.1725473e+00  -1.2535951e+01,4.5964720e-03,1.3025826e-03  -4.3892271e-01,1.2674668e+00,2.5691424e-01});
  }
  
  yasso15::yasso15(const std::array<double, 35> &atheta) : A{0.} {
    setTheta(atheta);
  }
  
  void yasso15::setTheta(const std::array<double, 35> &atheta) {
    for(int i=0; i<35; ++i) {theta[i] = atheta[i];}
    theta[31] = -std::fabs(theta[31]);
    theta[34] = -std::fabs(theta[34]);
    for(int i=0; i<4; ++i) {theta[i] = -std::fabs(theta[i]);}
    //aNeedToBeSet=true;
  }

  void yasso15::setClimSizeLeach(const double& avgT, const double& sumP, const double& ampT, const double& diam, const double& leach) {
    double tem{0.};
    double temN{0.};
    double temH{0.};
    const double m3 = sumP/1000.;
    {
      const double pi=M_PI;
      //temperature annual cycle approximation
      const double m0 = (1./sqrt(2.)-1.)/pi;
      const double m1 = 1./sqrt(2.)/pi;
      const double m2 = (1.-1./sqrt(2.))/pi;
      const double clim = 4. * ampT;
      double te[4];
      te[0] = avgT + clim * m0;
      te[1] = avgT - clim * m1;
      te[2] = avgT + clim * m2;
      te[3] = avgT + clim * m1;
      
      //Average temperature dependence
      double te2[4];
      for(int i=0; i<4; ++i) {te2[i] = te[i] * te[i];}
      for(int i=0; i<4; ++i) {
	tem += exp(theta[21]*te[i]+theta[22]*te2[i]);
	temN += exp(theta[23]*te[i]+theta[24]*te2[i]);
	temH += exp(theta[25]*te[i]+theta[26]*te2[i]);
      }
      tem /= 4.; temN /= 4.; temH /= 4.;
      
      //Precipitation dependence
      tem *= 1.-exp(theta[27] * m3);
      temN *= 1.-exp(theta[28] * m3);
      temH *= 1.-exp(theta[29] * m3);
    }
    
    //! check rare case where no decomposition happens for some compartments 
    //! (basically, if no rain)
    //gk: Changed that there is also a solution for spinnup
    if(tem < tol) {
      noDecomposition = true;
      return;
    }
    noDecomposition = false;
    
    //Size class dependence -- no effect if d == 0.0
    if(diam > 0.) { //gk: orig calculates also for negative diam
      const double size_dep = std::min(1., pow(1. + theta[32]*diam + theta[33]*diam*diam, theta[34]));
      //Calculating matrix a (will work ok despite the sign of alphas)
      for(int i=0; i<3; ++i) {A[i*6] = theta[i]*tem*size_dep;}
      A[3*6] = theta[3]*temN*size_dep;
    } else {
      for(int i=0; i<3; ++i) {A[i*6] = theta[i]*tem;}
      A[3*6] = theta[3]*temN;
    }
    A[24] = theta[31]*temH; //no size effect in humus
    const std::array<double,4> dAbs{std::fabs(A[0]), std::fabs(A[6]), std::fabs(A[2*6]), std::fabs(A[3*6])};
    for(int i=0, idx=3; i<4; ++i) {
      for(int j=0; j<4; ++j) {
	if(i!=j) {A[i*5+j] = theta[++idx] * dAbs[j];}
      }
      //gk: default 0 init:  A[i*5+4] = 0.; //no mass flows from H -> AWEN
    }
    //mass flows AWEN -> H (size effect is present here)
    for(int i=0; i<4; ++i) {A[20+i] = theta[30] * dAbs[i];}
    //Leaching (no leaching for humus) 
    if(leach < 0.) { //gk: orig calculates also for positive leach
      const double aux = leach * m3;
      for(int i=0; i<4; ++i) {A[6*i] += aux;}
    }
    //aNeedToBeSet=false;
    aNeedToBeDecomp=true;
    aNeedToBeExpo=true;
  }

  void yasso15::getSpin(const std::array<double, 5>& infall, std::array<double, 5>& result, const int years) {
    // if(aNeedToBeSet) {
    //   for(int i=0; i<5; ++i) {
    // 	result[i] = std::numeric_limits<double>::quiet_NaN();
    //   }
    // } else {
      if(noDecomposition) {
	for(int i=0; i<5; ++i) {
	  result[i] = years * infall[i];
	}
      } else {
	if(aNeedToBeDecomp) {
	  Crout<5>(A, Adecomp);
	  aNeedToBeDecomp=false;
	}
	solveCrout<5>(Adecomp, infall, result);
	for(int i=0; i<5; ++i) {result[i] *= -1.;}
      }
      //}
  }

  void yasso15::setTimespan(const double& atimespan) {
    timespan = atimespan;
    aNeedToBeExpo=true;
  }

  size_t yasso15::setTaylorTerms(const size_t& ataylorTerms) {
    taylorTerms = ataylorTerms+1;
    return(taylorTerms-1);
  }

  void yasso15::getNextTimestep(const std::array<double, 5>& init, const std::array<double, 5>& infall, std::array<double, 5>& result) {
    // if(aNeedToBeSet) {
    //   for(int i=0; i<5; ++i) {
    // 	result[i] = std::numeric_limits<double>::quiet_NaN();
    //   }
    // } else {
      if(noDecomposition) {
	for(int i=0; i<5; ++i) {result[i] += infall[i];}
	return;
      } else {
      //Solve the differential equation x'(t) = A(theta)*x(t) + b, x(0) = init
      //Solve DE in given time
      //Solve Matrix Differential equation Taylor x'(t) = Ax(t) + b
	std::array<double, 5> z1;
	MATMUL<5,5,1>(A,init,z1);
	for(int i=0; i<5; ++i) {z1[i] += infall[i];}
	if(aNeedToBeExpo) {
	  std::array<double, 5*5> At;
	  for(int i=0; i<5*5; ++i) {At[i] = A[i] * timespan;}
	  mexpAt.fill(0.);
	  matrixexp<5>(At, mexpAt, taylorTerms);
	  aNeedToBeExpo = false;
	}
	std::array<double, 5> z2;
	MATMUL<5,5,1>(mexpAt,z1,z2);
	for(int i=0; i<5; ++i) {z2[i] -= infall[i];}
	if(aNeedToBeDecomp) {
	  Crout<5>(A, Adecomp);
	  aNeedToBeDecomp = false;
	}
	solveCrout<5>(Adecomp, z2, result);
      }
      // }
  }

}
