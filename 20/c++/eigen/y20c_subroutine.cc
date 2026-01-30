#include "y20c_subroutine.h"
#include <Eigen/Dense>

namespace yasso {

using namespace Eigen;

typedef Matrix<double, 5, 5> Matrix5d;
typedef Matrix<double, 5, 1> Vector5d;
typedef Matrix<double, 35, 1> Vector35d;
typedef Matrix<double, 12, 1> Vector12d;

struct yasso20::Impl {
    Vector35d theta;
    Matrix5d A;
    Matrix5d mexpAt;

    Matrix5d exp5(const Matrix5d& mat, size_t q) {
        double p = mat.norm();
        double normiter = 2.0;
        int j = 2;
        
        while (p > normiter) {
            normiter *= 2.0;
            ++j;
        }

        Matrix5d C = mat / normiter; // Skalierung
        Matrix5d B = Matrix5d::Identity();
        B += C;
        
        Matrix5d D = C;
        
        // Taylor-Reihe bis q Terme
        for (size_t i = 2; i < q; ++i) {
            D = D * C;
            D /= static_cast<double>(i);
            B += D;
        }

        // Squaring-Phase
        for (int k = 1; k < j; ++k) {
            B = B * B;
        }
        
        return B;
    }

  Matrix5d exp5FIX(const Matrix5d& mat) {
        double s = mat.squaredNorm();
        double p = std::sqrt(s);
        double normiter = 2.0;
        int j = 2;
        while (p > normiter) {
            normiter *= 2.0;
            j += 1;
        }

        Matrix5d C = mat / normiter;
        Matrix5d res = Matrix5d::Identity() + C;
        Matrix5d D = C * C; res += D * 0.5;
        D *= C; res += D * 0.16666666666666666;
        D *= C; res += D * 0.041666666666666664;
        D *= C; res += D * 0.008333333333333333;
        D *= C; res += D * 0.001388888888888889;
        D *= C; res += D * 0.0001984126984126984;

        for (int k = 2; k <= j; ++k) res = res * res;
        return res;
    }


    void updateA(const std::array<double, 12>& avgT, double sumP, double diam, double leach) {
        double m3 = sumP / 1000.0;
        double tem_sum = 0, temN_sum = 0, temH_sum = 0;

        for (int i = 0; i < 12; ++i) {
            double t = avgT[i];
            double t2 = t * t;
            tem_sum  += std::exp(theta[21] * t + theta[22] * t2);
            temN_sum += std::exp(theta[23] * t + theta[24] * t2);
            temH_sum += std::exp(theta[25] * t + theta[26] * t2);
        }

        double tem  = (tem_sum / 12.0)  * (1.0 - std::exp(theta[27] * m3));
        double temN = (temN_sum / 12.0) * (1.0 - std::exp(theta[28] * m3));
        double temH = (temH_sum / 12.0) * (1.0 - std::exp(theta[29] * m3));

        double size_dep = (diam > 0.0) ? std::min(1.0, std::pow(1.0 + theta[32]*diam + theta[33]*diam*diam, theta[34])) : 1.0;

        A.setZero();
        double k = tem * size_dep;
        A(0,0) = theta[0] * k; A(1,1) = theta[1] * k; A(2,2) = theta[2] * k;
        A(3,3) = theta[3] * temN * size_dep;
        A(4,4) = theta[31] * temH;

        double dAbs[4] = { std::abs(A(0,0)), std::abs(A(1,1)), std::abs(A(2,2)), std::abs(A(3,3)) };
        int idx = 4;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (i != j) { A(i, j) = theta[idx] * dAbs[j]; idx++; }
            }
        }
        for (int j = 0; j < 4; ++j) A(4, j) = theta[30] * dAbs[j];

        if (leach < 0.0) {
            double lm3 = leach * m3;
            for (int j = 0; j < 4; ++j) A(j, j) += lm3;
        }
    }
};

yasso20::yasso20() : pimpl(std::make_unique<Impl>()) {}

yasso20::~yasso20() = default;

yasso20::yasso20(const std::array<double, 35> &t_arr) : yasso20() {
    setTheta(t_arr);
}

void yasso20::setTheta(const std::array<double, 35> &t_arr) {
    pimpl->theta = Vector35d::Map(t_arr.data());
    int tabs[] = {0, 1, 2, 3, 31, 34};
    for (int i : tabs) pimpl->theta[i] = -std::abs(pimpl->theta[i]);
}

void yasso20::setClimSizeLeach(const std::array<double, 12>& avgT, double sumP, double diam, double leach) {
    pimpl->updateA(avgT, sumP, diam, leach);
    noDecomposition = (pimpl->A.norm() < 1e-12);
    pimpl->mexpAt = pimpl->exp5(pimpl->A * timespan, taylorTerms);
}

void yasso20::getNextTimestep(const std::array<double, 5>& init, const std::array<double, 5>& infall, std::array<double, 5>& result) {
    Vector5d v_init = Vector5d::Map(init.data());
    Vector5d v_infall = Vector5d::Map(infall.data());
    
    Vector5d rhs = pimpl->mexpAt * (pimpl->A * v_init + v_infall) - v_infall;
    Vector5d res = pimpl->A.partialPivLu().solve(rhs);
    
    std::copy(res.data(), res.data() + 5, result.begin());
}

void yasso20::setTimespan(double t) { timespan = t; }

size_t yasso20::setTaylorTerms(size_t terms) {
    taylorTerms = terms;
    if (!noDecomposition) {
        pimpl->mexpAt = pimpl->exp5(pimpl->A * timespan, taylorTerms);
    }
    return taylorTerms;
}

void yasso20::getSpin(const std::array<double, 5>& infall, std::array<double, 5>& result, int years) {
    if (noDecomposition) {
        // Sonderfall: Kein Abbau -> Akkumulation über die Jahre
        for (int i = 0; i < 5; ++i) {
            result[i] = years * infall[i];
        }
    } else {
        // Stationärer Zustand (Steady State): A * x + infall = 0  =>  A * x = -infall
        Vector5d v_infall = Vector5d::Map(infall.data());
        Vector5d v_rhs = -v_infall; // Die rechte Seite ist -infall

        // Lösen des Systems A * x = v_rhs mittels LU-Zerlegung (entspricht funktional Crout)
        // PartialPivLU ist effizient und für 5x5 Matrizen sehr schnell.
        Vector5d x = pimpl->A.partialPivLu().solve(v_rhs);

        std::copy(x.data(), x.data() + 5, result.begin());
    }
}

} // namespace yasso
