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
  Matrix5d th;
  Matrix5d A;
  Matrix5d mexpAt;
  Eigen::PartialPivLU<Eigen::Matrix<double, 5, 5>> lu;

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
      double tem = 0, temN = 0, temH = 0, m3 = prec / 1000.;
        
        for (int i = 0; i < 12; ++i) {
            double t = temp[i];
            double t2 = t * t;
            tem  += std::exp(theta[21] * t + theta[22] * t2);
            temN += std::exp(theta[23] * t + theta[24] * t2);
            temH += std::exp(theta[25] * t + theta[26] * t2);
        }
        
        tem  = (tem / 12.0)  * (1.0 - std::exp(theta[27] * m3));
        temN = (temN / 12.0) * (1.0 - std::exp(theta[28] * m3));
        temH = (temH / 12.0) * (1.0 - std::exp(theta[29] * m3));

        double size_dep = std::min(1.0, std::pow(1.0 + theta(32)*d + theta(33)*d*d, theta(34)));

        if (tem <= 1e-12) [[unlikely]] {
	  A.setZero();
	  return;
	}

	Matrix<double, 5, 1> dA;
	dA.head<3>() = theta.head<3>() * (tem * size_dep);
	dA[3] = theta[3] * temN * size_dep;
	dA[4] = theta[31] * temH;
	//A = th; 
	//A.applyOnTheRight(dA.cwiseAbs().asDiagonal());
	A = th * dA.cwiseAbs().asDiagonal();
	A.diagonal() = dA;
 
        // Leaching
	if (leac < 0.) [[unlikely]]
	  A.diagonal().head<4>().array() += leac * m3;

	lu = A.partialPivLu();
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

    pimpl->th.diagonal().setConstant(-1.0);
    int i1[] = {1, 2, 3, 5, 7, 8, 10, 11, 13, 15, 16, 17};
    int i2[] = {7, 10, 13, 4, 11, 14, 5, 8, 15, 6, 9, 12};
    for (std::size_t i = 0; i < std::size(i1); ++i) {
      pimpl->th(i1[i]) = t_arr[i2[i]];
    }
    pimpl->th.row(4).head<4>().setConstant(t_arr[30]);
    pimpl->th.col(4).head<4>().setZero();
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
    //Vector5d xt = pimpl->A.partialPivLu().solve(z2);
    Vector5d xt = pimpl->lu.solve(z2);
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
