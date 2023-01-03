#include "y15cEigen_subroutine.h"

#include <unsupported/Eigen/MatrixFunctions>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
    aNeedToBeSet=true;
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
    for(int i=0; i<5; ++i) {
      for(int j=0; j<5; ++j) {
	Am(i,j) = A[i*5+j];
      }
    }
    aNeedToBeSet=false;
    aNeedToBeDecomp=true;
    aNeedToBeExpo=true;
  }

  void yasso15::getSpin(const Eigen::Matrix<double, 5, 1>& infall, Eigen::Matrix<double, 5, 1>& result) {
    if(aNeedToBeSet) {
      for(int i=0; i<5; ++i) {
	result[i] = std::numeric_limits<double>::quiet_NaN();
      }
    } else {
      if(noDecomposition) {
	for(int i=0; i<5; ++i) {
	  if(infall[i] > 0.) {result[i] = std::numeric_limits<double>::infinity();
	  } else {result[i] = 0.;}
	}
      } else {
	result = Am.householderQr().solve(infall) * -1.;
      }
    }
  }

  void yasso15::setTimespan(const double& atimespan) {
    timespan = atimespan;
    aNeedToBeExpo=true;
  }

  void yasso15::getNextTimestep(const Eigen::Matrix<double, 5, 1>& init, const Eigen::Matrix<double, 5, 1>& infall, Eigen::Matrix<double, 5, 1>& result) {
    if(aNeedToBeSet) {
      for(int i=0; i<5; ++i) {
	result[i] = std::numeric_limits<double>::quiet_NaN();
      }
    } else {
      if(noDecomposition) {
	for(int i=0; i<5; ++i) {result[i] += infall[i];}
	return;
      } else {
      //Solve the differential equation x'(t) = A(theta)*x(t) + b, x(0) = init
      //Solve DE in given time
      //Solve Matrix Differential equation Taylor x'(t) = Ax(t) + b
	z1 = Am * init + infall;
	z2 = (Am * timespan).exp() * z1 - infall;
	result = Am.householderQr().solve(z2);
      }
    }
  }

}
