#include <iostream>
#include <array>
#include <chrono>

#include "y20c_subroutine.h"

using namespace std;

int main() {
  { //Rewritten Yasso20
    array<double, 35> theta = {0.51,5.19,0.13,0.1,0.5,0.,1.,1.,0.99,0.,0.,0.,0.,0.,0.163,0.,-0.,0.,0.,0.,0.,0.158,-0.002,0.17,-0.005,0.067,-0.,-1.44,-2.0,-6.9,0.0042,0.0015,-2.55,1.24,0.25}; //Parameters
    double time {1.}; //Time to run
    array<double, 12> avgT {-2., 2., 5., 10., 15., 20., 22., 21., 16., 11., 6., 3.}; //Monthly mean temperatures [C]
    double sumP = 600.; //Precip annual summ [mm]
    //0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
    array<double, 5> init {0.,0.,0.,0.,0.}; //Initial State
    array<double, 5> infall {0.5,0.1,0.1,0.2,0.}; //Infall
    array<double, 5> result;
    double diam {2.}; //size [cm]
    double leach {0.}; //Leaching

    yasso::yasso20 yasso(theta);
    
    cout << "\nIn C++ rewritten Yasso20" << endl;
    yasso.setClimSizeLeach(avgT, sumP, diam, leach);
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
    array<double, 35> theta = {0.51,5.19,0.13,0.1,0.5,0.,1.,1.,0.99,0.,0.,0.,0.,0.,0.163,0.,-0.,0.,0.,0.,0.,0.158,-0.002,0.17,-0.005,0.067,-0.,-1.44,-2.0,-6.9,0.0042,0.0015,-2.55,1.24,0.25}; //Parameters
    double time {1.}; //Time to run
    array<double, 12> avgT {-2., 2., 5., 10., 15., 20., 22., 21., 16., 11., 6., 3.}; //Monthly mean temperatures [C]
    double sumP = 600.; //Precip annual summ [mm]
    //0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
    array<double, 5> init {0.,0.,0.,0.,0.}; //Initial State
    array<double, 5> infall {0.5,0.1,0.1,0.2,0.}; //Infall
    double diam {2.}; //size [cm]
    double leach {0.}; //Leaching
    array<double, 5> infall2;
    array<double, 12> avgT2;
    yasso::yasso20 yasso(theta);
    yasso.setTaylorTerms(6); //Default is 11

    auto start = std::chrono::high_resolution_clock::now();
    size_t n = 1000000;
    yasso.setTimespan(time);
    for(size_t i = 1; i <= n; ++i) {
      for(int j=0; j<5; ++j) {infall2[j] = infall[j] * double(i) / double(n);}
      for(int j=0; j<12; ++j) {avgT2[j] = avgT[j] + double(i) / double(n);}
      yasso.setClimSizeLeach(avgT2, sumP, diam, leach);
      yasso.getNextTimestep(init, infall2, init);
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    cout << "Time: " << elapsed.count() << endl;
    cout << "Result:";
    for(int j=0; j<5; ++j) {cout << " " << init[j];}
    cout << endl;
  }
  //Time: 0.549026
  //Result: 1.40537 0.149233 0.223336 3.02024 6.66991
  //Time: 0.38552 (with fastmath)
  
  return 0;
}
