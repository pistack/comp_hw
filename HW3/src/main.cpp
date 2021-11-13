/*!
 * @file main.cpp
 * @brief main program for homework3 of Computer1 class in Yonsei University
 * Interactively reads inital condition, order of basis function,
 * number of points to evaluate, number of iteration, step size, parameter lambda and
 * output file name then computes and saves solution.
 * @author pistack (Junho Lee)
 * @date 2021. 11. 12.
 */

#ifndef PATH_TYPE_FOURIER
#define PATH_TYPE_FOURIER // use default path type: fourier_path
#endif

#ifdef PATH_TYPE_BEZIER
#undef PATH_TYPE_FOURIER
#endif 

#ifndef PRECISION_LEVEL
#define PRECISION_LEVEL 1 // use default precision: double
#endif

#ifndef MONITOR
#define MONITOR 0 // default: does not monitor optimization process
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <chrono>
#include <string>
#include <iostream>
#include <fstream>

#ifdef PATH_TYPE_FOURIER
#include "libmcm/mcm_fourier.hpp"
#endif

#ifdef PATH_TYPE_BEZIER
#include "libmcm/mcm_bezier.hpp"
#endif

#if PRECISION_LEVEL == 0
    #define PRECISION float
    #define DIGITS 6
#elif PRECISION_LEVEL == 1
    #define PRECISION double
    #define DIGITS 14
#endif

using namespace std;

constexpr PRECISION pi = std::acos(PRECISION(-1)); // pi

/// @brief functor class for the kepler lagranian
/// @param T precision should be the one of
/// float, double, long double.

template<typename T>
class kepler_lag{

  public :

  kepler_lag() {}

  /// @brief evaluates kepler lagrangian at given time
  /// @param t time
  /// @param p path
  /// @param dp derivative of path
  /// @return lagrangian evaluated at given time
  T operator()(T t, 
  vector<T> p, vector<T> dp) const
  {
    if(p[0]<1e-2)
    return 0.5*(pow(dp[0], 2.0)+pow(p[0]*dp[1], 2.0))+100; // to avoid singularity and penalize bad guess
    return 0.5*(pow(dp[0], 2.0)+pow(p[0]*dp[1], 2.0))+
    1/p[0];
  }

};

int main(void)
{
  // execution time
  std::chrono::steady_clock::time_point start;
  std::chrono::steady_clock::time_point end;
  PRECISION rtol; // abs tol of action integral
  PRECISION e; // estimated error
  unsigned int num_eval; // number of points to eval
  unsigned int order; // order of basis function
  std::size_t max_iter; // number of iteration
  PRECISION zeta_min; // minimum value of zeta
  PRECISION t0; // initial time
  PRECISION max_step; // step size
  PRECISION lambda; // parameter to accept move
  PRECISION min_action; // minimum action value
  string filename; // file name to store results
  string filename_coeff; // file name to store coeffcients
  ofstream fout; // file output stream
  #ifdef MONITOR
  string filename_monitor; // file name to monitor optimization process
  #endif

  cout << "==========================================================" << endl;
  cout << "               hw3: main program for homework3            " << endl;
  cout << " initial time: ";
  cin >> t0;
  cout << " minimum value of zeta: ";
  cin >> zeta_min;
  cout << " relative tolerance of action: ";
  cin >> rtol;
  #ifdef PATH_TYPE_FOURIER
  cout << " number of sine and cosine function to use: ";
  #endif
  #ifdef PATH_TYPE_BEZIER
  cout << " order of bezier curve: ";
  #endif
  cin >> order;
  cout << " number of points to evaluate path: ";
  cin >> num_eval;
  cout << " size of step: ";
  cin >> max_step;
  cout << " value of paramter which controls acceptance of move: ";
  cin >> lambda;
  cout << " Number of iteration: ";
  cin >> max_iter;
  cout << " file name to store result: ";
  cin >> filename;
  cout << " file name to store coefficient: ";
  cin >> filename_coeff;
  #if MONITOR == 1
  cout << " file name to monitor optimization process: ";
  cin >> filename_monitor;
  #endif
  cout << " Now starts calculation" << endl;
  start = std::chrono::steady_clock::now();
  // initial condition
  PRECISION zeta_max = zeta_min/(2*zeta_min-1);
  PRECISION a = (zeta_min+zeta_max)/2;
  PRECISION tmax = pi*pow(a, 1.5);
  vector<PRECISION> p0 = {zeta_min, 0};
  vector<PRECISION> p1 = {zeta_max, pi};
  #ifdef PATH_TYPE_FOURIER
  PRECISION add_setup = 2*tmax;
  libmcm::mcm_fourier<PRECISION, 
  kepler_lag<PRECISION>> kepler(PRECISION(0), 
  tmax, p0, p1, rtol, order, add_setup);
  #endif

  #ifdef PATH_TYPE_BEZIER
  libmcm::mcm_bezier<PRECISION, 
  kepler_lag<PRECISION>> kepler(PRECISION(0), 
  tmax, p0, p1, rtol, order);
  #endif

  kepler.set_init_guess();
  // variable to get optimization state
  int num_move; // number of actual moves
  PRECISION accept_ratio; // acceptance ratio

  #if MONITOR == 1
  tie(num_move, accept_ratio) = \
  kepler.optimize(max_iter, max_step, lambda, 
  filename_monitor);
  #else
  tie(num_move, accept_ratio) = \
  kepler.optimize(max_iter, max_step, lambda);
  #endif

  // Now fill time
  vector<PRECISION> t(num_eval+1, 0); // one more point need for end point
  PRECISION grid_space = tmax/PRECISION(num_eval);
  for(int i=1; i<=num_eval; i++)
  t[i] = t[i-1]+grid_space;
  t[num_eval] = tmax;

  // store result
  vector<vector<PRECISION>> result(2, vector<PRECISION>(num_eval+1, 0));
  result = kepler.min_eval(t);
  min_action = kepler.get_min_action(e);

  // store coefficients
  // fourier
  #ifdef PATH_TYPE_FOURIER
  vector<PRECISION> adder;
  vector<PRECISION> scaler;
  vector<vector<PRECISION>> coeff;

  tie(adder, scaler, coeff) = kepler.get_min_coeff();

  unsigned int dim1 = adder.size();
  unsigned int dim2 = coeff[0].size();
  #endif

  // bezier
  #ifdef PATH_TYPE_BEZIER
  vector<PRECISION> scaler;
  vector<vector<PRECISION>> coeff;

  tie(scaler, coeff) = kepler.get_min_coeff();

  unsigned int dim1 = scaler.size();
  unsigned int dim2 = coeff[0].size();
  #endif

  // move t by t0
  transform(t.begin(), t.end(), t.begin(),
  [t0](PRECISION &x){return x += t0;});
  end = std::chrono::steady_clock::now();
  cout << " Calcuation is finished" << endl;
  cout.unsetf(ios::floatfield); // initialize floatfield
  cout << fixed;
  cout.precision(2); // print significant digits
  cout << " Elapsed time: " << \
  std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000000.0
  << " second" << endl;
  cout << "======================result==============================" << endl;
  cout.unsetf(ios::floatfield); // initialize floatfield
  cout.precision(DIGITS); // print significant digits
  cout << " Minimum action is " << min_action << endl;
  cout << " Estimated error of action integration is " << e <<endl;
  cout << " Number of Actual move is " << num_move << endl;
  cout << " Acceptance ratio: " << accept_ratio << endl;
  // store results to file
  fout.open(filename);
  fout << '#' << 't' << '\t' << "zeta" << '\t' << "theta" << endl;
  fout.unsetf(ios::floatfield); // initialize floatfield
  fout.precision(DIGITS); // print significant digits
  for(unsigned int i=0; i < num_eval+1; i++)
    fout << t[i] << '\t' << result[0][i] << '\t' << result[1][i] << endl;
  fout.close();
  fout.open(filename_coeff);
  fout << '#' << '\t' << "zeta" << '\t' << "theta" << endl;
  #ifdef PATH_TYPE_FOURIER
  fout << "ADDER" << endl;
  fout << adder[0] << '\t' << adder[1] << endl;
  #endif
  fout << "SCALER" << endl;
  fout << scaler[0] << '\t' << scaler[1] << endl;
  fout << "COEFF" << endl;
  for(unsigned int i=0; i<dim2; ++i)
  fout << coeff[0][i] << '\t' << coeff[1][i] << endl;
  fout.close(); 
  cout << " Save result to " << filename << endl;
  cout << " Save coeffcients to " << filename_coeff << endl;
  cout << " Teriminates program, good bye :) " << endl;
  cout << "==========================================================" << endl;

  return 0;
}
      
      

  
						  
