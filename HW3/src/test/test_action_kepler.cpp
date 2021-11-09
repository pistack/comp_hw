/*!
 * @file test_action_kepler.cpp
 * @brief test action::eval() routine with kepler action
 * @author pistack (Junho Lee)
 * @date 2021. 11. 9.
 */

#include <cmath>
#include <chrono>
#include <iostream>
#include "fourier_path.hpp"
#include "action.hpp"

#if PRECISION_LEVEL == 0
    #define PRECISION float
    #define DIGITS 6
#elif PRECISION_LEVEL == 1
    #define PRECISION double
    #define DIGITS 14
#endif

using namespace libpath;
using namespace std;

#ifndef DOXYGEN_SKIP
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
#endif /* DOXYGEN_SKIP */

int main(void)
{

  cout << "==========================================================" << endl;
  cout << "               Test action::eval() routine                " << endl;
  cout << " Test 1. n_f: 3 with kepler action                        " << endl;
  cout << " Test 2. n_f: 4 with kepler action                        " << endl;

  // measure execution time of action::eval() method
  std::chrono::steady_clock::time_point start;
  std::chrono::steady_clock::time_point end;
  // estimated error
  PRECISION e;

  // initial condition
  PRECISION zeta_min = 0.9;
  PRECISION zeta_max = zeta_min/(2*zeta_min-1);
  PRECISION a = 0.5*(zeta_min+zeta_max);
  PRECISION tmax = pi<PRECISION>*pow(a, 1.5);
  PRECISION period = 2*tmax;
  #if PRECISION_LEVEL == 0
  vector<PRECISION> tol = {1.0, 1e-2, 1e-4};
  #elif PRECISION_LEVEL == 1
  vector<PRECISION> tol = {1.0, 1e-2, 1e-4, 1e-6, 1e-8};
  #endif
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
    fourier<PRECISION>(3, period, c1[0])),
    fourier_path<PRECISION>(0.0, tmax, 0.0, pi<PRECISION>,
    fourier<PRECISION>(3, period, c1[1])),
  };
  vector<fourier_path<PRECISION>> path2 {
    fourier_path<PRECISION>(0.0, tmax, zeta_min, zeta_max,
    fourier<PRECISION>(4, period, c2[0])),
    fourier_path<PRECISION>(0.0, tmax, 0.0, pi<PRECISION>,
    fourier<PRECISION>(4, period, c2[1])),
  };
  action<PRECISION, fourier_path<PRECISION>, kepler_lag<PRECISION>> tst1(path1);
  action<PRECISION, fourier_path<PRECISION>, kepler_lag<PRECISION>> tst2(path2);

  cout.unsetf(ios::floatfield); // initialize floatfield
  cout.precision(DIGITS); // print significant digits
  for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst1.update(*it);
    cout << " Test 1. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: gauss-kronrod quadrature" << endl;
    cout << "Integration value: " << tst1.eval(e) << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    cout << "Integration Method: Tanh-Sinh quadrature" << endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration value: " << tst1.eval(1, 7, e) << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst2.update(*it);
    cout << " Test 2. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst2.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: gauss-kronrod quadrature" << endl;
    cout << "Integration value: " << tst2.eval(e) << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    cout << "Integration Method: Tanh-Sinh quadrature" << endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst2.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration value: " << tst2.eval(1, 7, e) << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  cout << "Test finished!" << endl;
  cout << "==========================================================" << endl;

  return 0;
}
      
      

  
						  
