#include <iostream>
#include <array>
#include <chrono>

using namespace std;

extern"C" {
  void __yasso20_MOD_mod5c20(float *theta, float *time, float *temp, float *prec, float *init, float *b, float *d, float *leach, float *xt, int *steadystate_pred);
}

int main() {
  { //Yasso20 Fortran
    float theta[35] = {0.51,5.19,0.13,0.1,0.5,0.,1.,1.,0.99,0.,0.,0.,0.,0.,0.163,0.,-0.,0.,0.,0.,0.,0.158,-0.002,0.17,-0.005,0.067,-0.,-1.44,-2.0,-6.9,0.0042,0.0015,-2.55,1.24,0.25}; //Parameters
    float time {1.}; //Time to run
    float temp[12] {-2., 2., 5., 10., 15., 20., 22., 21., 16., 11., 6., 3.}; //Monthly mean temperatures [C]
    float precip {600.}; //Precipitation annual summ [mm]
    //0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
    float init[5] {0.,0.,0.,0.,0.}; //Initial State
    float b[5] {0.5,0.1,0.1,0.2,0.}; //Infall
    float xt[5] {0.,0.,0.,0.,0.}; //Result
    float d {2.}; //size [cm]
    float leach {0.}; //Leaching
    int steadystate_pred {0}; //set to true if ignore 'time' and compute solution in steady-state conditions (which sould give equal solution as if time is set large enough)
    
    cout << "Fortran Yasso20" << endl;
    for(int year=0; year<10; ++year) {
      __yasso20_MOD_mod5c20(theta, &time, temp, &precip, init, b, &d, &leach, xt, &steadystate_pred);
      cout << year;
      for(int i=0; i<5; ++i) {cout << " " << xt[i];}
      cout << endl;
      for(int i=0; i<5; ++i) {init[i] = xt[i];}
    }
    steadystate_pred = 1;
    __yasso20_MOD_mod5c20(theta, &time, temp, &precip, init, b, &d, &leach, xt, &steadystate_pred);
    cout << "*";
    for(int i=0; i<5; ++i) {cout << " " << xt[i];}
    cout << endl;
  }

  {
    float theta[35] = {0.51,5.19,0.13,0.1,0.5,0.,1.,1.,0.99,0.,0.,0.,0.,0.,0.163,0.,-0.,0.,0.,0.,0.,0.158,-0.002,0.17,-0.005,0.067,-0.,-1.44,-2.0,-6.9,0.0042,0.0015,-2.55,1.24,0.25}; //Parameters
    float time {1.}; //Time to run
    float temp[12] {-2., 2., 5., 10., 15., 20., 22., 21., 16., 11., 6., 3.}; //Monthly mean temperatures [C]
    float precip {600.}; //Precipitation annual summ [mm]
    //0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
    float init[5] {0.,0.,0.,0.,0.}; //Initial State
    float b[5] {0.5,0.1,0.1,0.2,0.}; //Infall
    float xt[5] {0.,0.,0.,0.,0.}; //Result
    float d {2.}; //size [cm]
    float leach {0.}; //Leaching
    float b2[5];
    float temp2[12];
    int steadystate_pred {0};

    auto start = std::chrono::high_resolution_clock::now();
    size_t n = 1000000;
    for(size_t i = 1; i <= n; ++i) {
      for(int j=0; j<5; ++j) {b2[j] = b[j] * double(i) / double(n);}
      for(int j=0; j<12; ++j) {temp2[j] = temp[j] + double(i) / double(n);}
      __yasso20_MOD_mod5c20(theta, &time, temp2, &precip, init, b2, &d, &leach, xt, &steadystate_pred);
      for(int i=0; i<5; ++i) {init[i] = xt[i];}
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    cout << "Time: " << elapsed.count() << endl;
    cout << "Result:";
    for(int j=0; j<5; ++j) {cout << " " << xt[j];}
    cout << endl;
  }
  //Time: 0.608852
  //Result: 1.40537 0.149233 0.223336 3.02024 6.66991

  return 0;
}
