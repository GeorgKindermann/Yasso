#include <iostream>
#include <array>

using namespace std;

extern"C" {
  void __yasso07_MOD_mod5c(float *theta, float *time, float *climate, float *init, float *b, float *d, float *leach, float *xt);
}

int main() {
  { //Yasso07 Fortran
    float theta[27] = {-0.49,-4.9,-0.24,-0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,0.091,-0.00021,-1.8,-0.15,-0.02,-0.92,-0.44,1.3,-0.26,-0.0013,0.0046};
    float time {1.}; //Time to run
    float climate[3] {10., 600., 12.}; //Temp annual average [C], precip annual summ [mm], amplitude (max. difference of month averages / 2) [C]
    //0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
    float init[5] {0.,0.,0.,0.,0.}; //Initial State
    float b[5] {0.5,0.1,0.1,0.2,0.}; //Infall
    float xt[5] {0.,0.,0.,0.,0.}; //Result
    float d {2.}; //size [cm]
    float leach {0.}; //Leaching

    cout << "Yasso07 Fortran" << endl;
    for(int year=0; year<10; ++year) {
      __yasso07_MOD_mod5c(theta, &time, climate, init, b, &d, &leach, xt);
      cout << year;
      for(int i=0; i<5; ++i) {cout << " " << xt[i];}
      cout << endl;
      for(int i=0; i<5; ++i) {init[i] = xt[i];}
    }
    time = 9999.;
    __yasso07_MOD_mod5c(theta, &time, climate, init, b, &d, &leach, xt);
    cout << "*";
    for(int i=0; i<5; ++i) {cout << " " << xt[i];}
    cout << endl;

  }

  return 0;
}
