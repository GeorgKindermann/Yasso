#include "y20c_subroutine.h"
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>

namespace yasso {

using namespace Eigen;

typedef Matrix<double, 5, 5> Matrix5d;
typedef Matrix<double, 5, 1> Vector5d;
typedef Matrix<double, 35, 1> Vector35d;

struct yasso20::Impl {
    Vector35d theta;
    Matrix5d A;
    Matrix5d mexpAt;

    Matrix5d matrixexp(const Matrix5d& At, size_t q) {
        double p = At.norm(); 
        double normiter = 2.0;
        int j = 1;
        
        while (p >= normiter) {
            normiter *= 2.0;
            j++;
        }

        Matrix5d C = At / normiter;
        Matrix5d B = Matrix5d::Identity();
        B += C;
        
        Matrix5d D = C;
        for (size_t i = 2; i < q; ++i) {
          D = (C * D) / static_cast<double>(i);
          B += D;
        }

        for (int k = 0; k < j; ++k) {
          B = B * B;
        }
        return B;
    }

    void updateA(double prec, double d, double leac, const std::array<double, 12>& temp) {
        double tem = 0, temN = 0, temH = 0;
        
        for (int i = 0; i < 12; ++i) {
            double t = temp[i];
            double t2 = t * t;
            tem  += std::exp(theta[21] * t + theta[22] * t2);
            temN += std::exp(theta[23] * t + theta[24] * t2);
            temH += std::exp(theta[25] * t + theta[26] * t2);
        }
        
        tem  = (tem / 12.0)  * (1.0 - std::exp(theta[27] * prec / 1000.0));
        temN = (temN / 12.0) * (1.0 - std::exp(theta[28] * prec / 1000.0));
        temH = (temH / 12.0) * (1.0 - std::exp(theta[29] * prec / 1000.0));

        double size_dep = std::min(1.0, std::pow(1.0 + theta(32)*d + theta(33)*d*d, theta(34)));


        A.setZero();
        if (tem <= 1e-12) return;

        // Diagonalelemente
        for(int i=0; i<3; ++i) A(i,i) = -std::abs(theta[i]) * tem * size_dep;
        A(3,3) = -std::abs(theta[3]) * temN * size_dep;

        // Massenflüsse
        A(0,1) = theta[4] * std::abs(A(1,1)); A(0,2) = theta[5] * std::abs(A(2,2)); A(0,3) = theta[6] * std::abs(A(3,3));
        A(1,0) = theta[7] * std::abs(A(0,0)); A(1,2) = theta[8] * std::abs(A(2,2)); A(1,3) = theta[9] * std::abs(A(3,3));
        A(2,0) = theta[10]* std::abs(A(0,0)); A(2,1) = theta[11]* std::abs(A(1,1)); A(2,3) = theta[12]* std::abs(A(3,3));
        A(3,0) = theta[13]* std::abs(A(0,0)); A(3,1) = theta[14]* std::abs(A(1,1)); A(3,2) = theta[15]* std::abs(A(2,2));
        
        // Humus
        A(4,4) = -std::abs(theta[31]) * temH;
        for(int i=0; i<4; ++i) A(4,i) = theta(30) * std::abs(A(i,i));

        // Leaching
        for(int i=0; i<4; ++i) A(i,i) += leac * prec / 1000.0;

    }
};

yasso20::yasso20() : pimpl(std::make_unique<Impl>()) {}

yasso20::yasso20(const std::array<double, 35> &t_arr) : yasso20() {
    setTheta(t_arr);
}

yasso20::~yasso20() = default;

void yasso20::setTheta(const std::array<double, 35> &t_arr) {
    pimpl->theta = Vector35d::Map(t_arr.data());
    
    // Die spezifischen Abbau-Raten (alpha) negativ forcieren
    int alpha_indices[] = {0, 1, 2, 3, 31, 34};
    for (int i : alpha_indices) {
        pimpl->theta[i] = -std::abs(pimpl->theta[i]);
    }
}

void yasso20::setClimSizeLeach(const std::array<double, 12>& avgT, double sumP, double diam, double leach) {
    pimpl->updateA(sumP, diam, leach, avgT);
    noDecomposition = (pimpl->A.norm() < 1e-12);
    if (!noDecomposition) {
      pimpl->mexpAt = pimpl->matrixexp(pimpl->A * timespan, taylorTerms);
    }
}

void yasso20::getNextTimestep(const std::array<double, 5>& init, const std::array<double, 5>& infall, std::array<double, 5>& result) {
    if (noDecomposition) {
        for(int i=0; i<5; ++i) result[i] = init[i] + infall[i] * timespan;
        return;
    }
    Vector5d v_init = Vector5d::Map(init.data());
    Vector5d v_b = Vector5d::Map(infall.data());
    Vector5d z1 = pimpl->A * v_init + v_b;
    Vector5d z2 = (pimpl->mexpAt * z1) - v_b;
    Vector5d xt = pimpl->A.partialPivLu().solve(z2);
    std::copy(xt.data(), xt.data() + 5, result.begin());
}

void yasso20::getSpin(const std::array<double, 5>& infall, std::array<double, 5>& result, int years) {
    if (noDecomposition) {
        for(int i=0; i<5; ++i) result[i] = (double)years * infall[i];
    } else {
        Vector5d v_b = Vector5d::Map(infall.data());
        Vector5d xt = (-pimpl->A).partialPivLu().solve(v_b);
        std::copy(xt.data(), xt.data() + 5, result.begin());
    }
}

size_t yasso20::setTaylorTerms(size_t terms) {
  taylorTerms = terms;
  if (!noDecomposition) {
    pimpl->mexpAt = pimpl->matrixexp(pimpl->A * timespan, taylorTerms);
  }
  return taylorTerms;
}

void yasso20::setTimespan(double t) { timespan = t; }

} // namespace yasso
