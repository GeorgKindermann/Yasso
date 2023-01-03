#define VARIATECLIMA 1
#define VARIATEINFALL 1
//#define SPINNUP
//VARIATECLIMA 0..no, 1..yes
//VARIATEINFALL 0..no, 1..yes

#include <iostream>
#include <array>
#include <chrono>
#include <armadillo>

#include "y15cArma_subroutine.h"

using namespace std;

int main() {
  size_t n = 1000000;
  cout << "Number of loops: " << n << endl;
  cout << "Variet Clima: " << VARIATECLIMA << endl;
  cout << "Variet Infall: " << VARIATEINFALL << endl;
  
  cout << "CPP Version" << endl;
  array<double, 35> theta = {0.49,4.9,0.24,0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,-0.15,-0.02,-0.92,-0.0004,-0.00017,0.091,-0.00021,0.049,-7.90E-05,0.035,-0.00021,-1.8,-1.2,-13,0.0046,0.0013,-0.44,1.3,0.26};
  double time {1.};
  double avgT = 10.;
  double sumP = 600.;
  double ampT = 12.;
  arma::vec5 init {0.,0.,0.,0.,0.};
  arma::vec5 infall {0.5,0.1,0.1,0.2,0.};
  arma::vec5 result;
  double diam {2.};
  double leach {0.};
  arma::vec5 infall2 {0.5,0.1,0.1,0.2,0.};
  double avgT2 = 10.;
  yasso::yasso15 yasso(theta);
  yasso.setTimespan(time);
  yasso.setClimSizeLeach(avgT2, sumP, ampT, diam, leach);

#ifdef SPINNUP
  cout << "\nSpinnup" << endl;
#else
  cout << "\nTimeserie" << endl;
#endif
  array<long double, 5> sum {0.,0.,0.,0.,0.};
  auto start = std::chrono::high_resolution_clock::now();
  for(size_t i = 0; i < n; ++i) {
#if VARIATEINFALL == 1
    for(int j=0; j<5; ++j) {infall2[j] = infall[j] * double(i) / double(n);}
#endif
#if VARIATECLIMA == 1
    avgT2 = avgT + double(i) / double(n);
    yasso.setClimSizeLeach(avgT2, sumP, ampT, diam, leach);
#endif
#ifdef SPINNUP
    yasso.getSpin(infall2, result);
    for(int j=0; j<5; ++j) {sum[j] += result[j];}
#else
    yasso.getNextTimestep(init, infall2, init);
    for(int j=0; j<5; ++j) {sum[j] += init[j];}
#endif
  }
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  cout << "Time: " << elapsed.count() << endl;
  cout << "Result:";
  for(int j=0; j<5; ++j) {cout << " " << sum[j] / double(n);}
  cout << endl;

  return(0);
}
