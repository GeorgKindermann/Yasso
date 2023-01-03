#include <iostream>
#include <array>
#include <chrono>

#include "y15c_subroutine.h"

using namespace std;

extern"C" {
  void __yasso_MOD_mod5c(float *theta, float *time, float *climate, float *init, float *b, float *d, float *leach, float *xt, int *steadystate_pred);
}

int main() {
  { //Fortran Yasso15
    float theta[35] = {0.49,4.9,0.24,0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,-0.15,-0.02,-0.92,-0.0004,-0.00017,0.091,-0.00021,0.049,-7.90E-05,0.035,-0.00021,-1.8,-1.2,-13,0.0046,0.0013,-0.44,1.3,0.26}; //Parameters
    float time {1.}; //Time to run
    float climate[3] {10., 600., 12.}; //Temp annual average [C], precip annual summ [mm], amplitude (max. difference of month averages / 2) [C]
    //0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
    float init[5] {0.,0.,0.,0.,0.}; //Initial State
    float b[5] {0.5,0.1,0.1,0.2,0.}; //Infall
    float xt[5] {0.,0.,0.,0.,0.}; //Result
    float d {2.}; //size [cm]
    float leach {0.}; //Leaching
    int steadystate_pred {0}; //set to true if ignore 'time' and compute solution in steady-state conditions (which sould give equal solution as if time is set large enough)
    
    cout << "Fortran Yasso15" << endl;
    for(int year=0; year<10; ++year) {
      __yasso_MOD_mod5c(theta, &time, climate, init, b, &d, &leach, xt, &steadystate_pred);
      cout << year;
      for(int i=0; i<5; ++i) {cout << " " << xt[i];}
      cout << endl;
      for(int i=0; i<5; ++i) {init[i] = xt[i];}
    }
    steadystate_pred = 1;
    __yasso_MOD_mod5c(theta, &time, climate, init, b, &d, &leach, xt, &steadystate_pred);
    cout << "*";
    for(int i=0; i<5; ++i) {cout << " " << xt[i];}
    cout << endl;
  }

  { //Rewritten Yasso15
    array<double, 35> theta = {0.49,4.9,0.24,0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,-0.15,-0.02,-0.92,-0.0004,-0.00017,0.091,-0.00021,0.049,-7.90E-05,0.035,-0.00021,-1.8,-1.2,-13,0.0046,0.0013,-0.44,1.3,0.26}; //Parameters
    double time {1.}; //Time to run
    double avgT = 10.; //Temp annual average [C]
    double sumP = 600.; //Precip annual summ [mm]
    double ampT = 12.; //Amplitude (max. difference of month averages / 2) [C]
    //0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
    array<double, 5> init {0.,0.,0.,0.,0.}; //Initial State
    array<double, 5> infall {0.5,0.1,0.1,0.2,0.}; //Infall
    array<double, 5> result;
    double diam {2.}; //size [cm]
    double leach {0.}; //Leaching

    yasso::yasso15 yasso(theta);
    
    cout << "\nIn C++ rewritten Yasso15" << endl;
    yasso.setClimSizeLeach(avgT, sumP, ampT, diam, leach);
    yasso.setTimespan(time);
    for(int year=0; year<10; ++year) {
      yasso.getNextTimestep(init, infall, init); //Write result in init
      cout << year;
      for(int i=0; i<5; ++i) {cout << " " << init[i];}
      cout << endl;
    }
    yasso.getSpin(infall, result);
    cout << "*";
    for(int i=0; i<5; ++i) {cout << " " << result[i];}
    cout << endl;
  }

  {
    array<double, 35> theta = {0.49,4.9,0.24,0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,-0.15,-0.02,-0.92,-0.0004,-0.00017,0.091,-0.00021,0.049,-7.90E-05,0.035,-0.00021,-1.8,-1.2,-13,0.0046,0.0013,-0.44,1.3,0.26}; //Parameters
    double time {1.}; //Time to run
    double avgT = 10.; //Temp annual average [C]
    double sumP = 600.; //Precip annual summ [mm]
    double ampT = 12.; //Amplitude (max. difference of month averages / 2) [C]
    //0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
    array<double, 5> init {0.,0.,0.,0.,0.}; //Initial State
    array<double, 5> infall {0.5,0.1,0.1,0.2,0.}; //Infall
    double diam {2.}; //size [cm]
    double leach {0.}; //Leaching
    array<double, 5> infall2 {0.5,0.1,0.1,0.2,0.};
    double avgT2 = 10.;
    yasso::yasso15 yasso(theta);
    yasso.setTaylorTerms(6); //Default is 11
 
    auto start = std::chrono::high_resolution_clock::now();
    size_t n = 1000000;
    yasso.setTimespan(time);
    for(size_t i = 1; i <= n; ++i) {
      for(int j=0; j<5; ++j) {infall2[j] = infall[j] * double(i) / double(n);}
      avgT2 = avgT + double(i) / double(n);
      yasso.setClimSizeLeach(avgT2, sumP, ampT, diam, leach);
      yasso.getNextTimestep(init, infall2, init);
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    cout << "Time: " << elapsed.count() << endl;
    cout << "Result:";
    for(int j=0; j<5; ++j) {cout << " " << init[j];}
    cout << endl;
  }
  //Time: 0.394529
  //Result: 2.61321 0.275452 0.392322 8.35294 10.5601

  return(0);
}
