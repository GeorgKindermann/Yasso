#ifndef YASSO_H_ARMA
#define YASSO_H_ARMA

//Rewritten Yasso15 core code based on Fortran code from Mikko Tuomi and Marko Järvenpää

#include <array>
#include <cmath>
#include <limits>
#include <armadillo>

namespace yasso {

  class yasso15 {
  public:
    yasso15();
    //Parameters - Theta
 // 1-16 matrix A entries: 4*alpha, 12*p
 // 17-21 Leaching parameters: w1,...,w5 IGNORED IN THIS FUNCTION
 // 22-23 Temperature-dependence parameters for AWE fractions: beta_1, beta_2
 // 24-25 Temperature-dependence parameters for N fraction: beta_N1, beta_N2
 // 26-27 Temperature-dependence parameters for H fraction: beta_H1, beta_H2
 // 28-30 Precipitation-dependence parameters for AWE, N and H fraction: gamma, gamma_N, gamma_H
 // 31-32 Humus decomposition parameters: p_H, alpha_H (Note the order!)
 // 33-35 Woody parameters: theta_1, theta_2, r
    yasso15(const std::array<double, 35> &theta);
    void setTheta(const std::array<double, 35> &theta);
 //avgT .. Temp annual average [C]
 //sumP .. Precip annual summ [mm]
 //ampT .. Amplitude (max. difference of month averages / 2) [C]
 //diam .. size [cm]
 //leach .. Leaching - Value below 0
    void setClimSizeLeach(const double& avgT, const double& sumP, const double& ampT, const double& diam, const double& leach);
 //timespann .. Time to run the model for one time step
    void setTimespan(const double& timespan);
    bool isThereDecomposition() {return(!noDecomposition);}
 //0..Ethanol sol.(waxes), 1..Water sol.(sugar), 2..Acid soluble (Celluloses), 3..Insoluble (lignin), 4..Humus : in init, infall an result
    void getSpin(const arma::vec5& infall, arma::vec5& result);
    void getNextTimestep(const arma::vec5& init, const arma::vec5& infall, arma::vec5& result);
  private:
    double timespan = 1.;
    std::array<double, 35> theta;
    std::array<double, 5*5> A;
    arma::mat55 Am;
    double tol = 1.E-12;
    arma::vec5 z1, z2;
    bool aNeedToBeSet=true;
    bool aNeedToBeDecomp=true;
    bool aNeedToBeExpo=true;
    bool noDecomposition = true;
  };
  
}

#endif
