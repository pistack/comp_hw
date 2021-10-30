/*!
 * @file test_action_simple.cpp
 * @brief test action::eval() routine with kepler action
 * @author pistack (Junho Lee)
 * @date 2021. 10. 30.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include "action.hpp"

#if PRECISION_LEVEL == 0
    #define PRECISION float
    #define DIGITS 6
#elif PRECISION_LEVEL == 1
    #define PRECISION double
    #define DIGITS 14
#endif

const PRECISION pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062; // define pi

using namespace std;

template<typename T>
class kepler_lag{

  public :

  kepler_lag() {}

  T operator()(T t, 
  vector<T> p, vector<T> dp) const
  {
    return 0.5*(pow(dp[0], 2.0)+pow(p[0]*dp[1], 2.0))+
    1/abs(p[0]);
  }

};

int main(void)
{

  cout << "==========================================================" << endl;
  cout << "               Test eval action routine                   " << endl;
  cout << " Test 1. single precision systematic approach n_f: 3      " << endl;
  cout << " Test 2. double precision n_f: 4 with atol: 1e-4          " << endl;


  // initial condition
  PRECISION zeta_min = 0.9;
  PRECISION zeta_max = zeta_min/(2*zeta_min-1);
  PRECISION a = 0.5*(zeta_min+zeta_max);
  PRECISION tmax = pi*pow(a, 1.5);
  PRECISION period = 2*tmax;
  vector<PRECISION> tol = {1e-4, 1e-6, 1e-8};
  vector<vector<PRECISION>> c1 {
    {-0.791412, -0.62261, -0.859226, -0.856224, 0.143711, -0.542851},
    {-0.216721, 0.548709, -0.198616, -0.0924677, 0.0491277, -0.0405731}
    };
  vector<vector<PRECISION>> c2 {
    {-0.47023040093914, -0.89070790762186, -0.11280783300074, -0.41071011096262,
    0.092578934740922, -0.055782211762309, 0.00061565328068328, -0.047040461332351},
    {-0.19905402253529, 0.70514094773255, -0.09667403548441, -0.14526040902903,
    0.12644790407768, 0.072545362434883, -0.045330151116832, 0.041103268443113}
    };
  vector<fourier_path<PRECISION>> path1 {
    fourier_path<PRECISION>(0.0, tmax, zeta_min, zeta_max,
    fourier<PRECISION>(3, 2*tmax, c1[0])),
    fourier_path<PRECISION>(0.0, tmax, 0.0, pi,
    fourier<PRECISION>(3, 2*tmax, c1[1])),
  };
  vector<fourier_path<PRECISION>> path2 {
    fourier_path<PRECISION>(0.0, tmax, zeta_min, zeta_max,
    fourier<PRECISION>(4, 2*tmax, c2[0])),
    fourier_path<PRECISION>(0.0, tmax, 0.0, pi,
    fourier<PRECISION>(4, 2*tmax, c2[1])),
  };
  action<PRECISION, kepler_lag<PRECISION>> tst1(path1);
  action<PRECISION, kepler_lag<PRECISION>> tst2(path2);

  cout.unsetf(ios::floatfield); // initialize floatfield
  cout.precision(DIGITS); // print significant digits
  for(int i=0; i<3; i++)
  {
    tst1.update(tol[i]);
    cout << " Test 1. atol: " << tol[i] <<  endl;
    cout << "Integration value: " << tst1.eval() << endl;
  }
  for(int i=0; i<3; i++)
  {
    tst2.update(tol[i]);
    cout << " Test 2. atol: " << tol[i] <<  endl;
    cout << "Integration value: " << tst2.eval() << endl;
  }
  cout << "Test finished!" << endl;
  cout << "==========================================================" << endl;

  return 0;
}
      
      

  
						  
