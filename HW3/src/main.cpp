/*!
 * @file main.cpp
 * @brief main program for homework3 of Computer1 class in Yonsei University
 * Interactively reads inital condition, number of sine function used for guess,
 * number of points to evaluate, number of interation, step size and
 * output file name then computes and saves solution.
 * @author pistack (Junho Lee)
 * @date 2021. 11. 04.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include "mcm.hpp"

#if PRECISION_LEVEL == 0
    #define PRECISION float
    #define DIGITS 6
#elif PRECISION_LEVEL == 1
    #define PRECISION double
    #define DIGITS 14
#endif

using namespace std;

/// @brief functor class for the kepler lagranian
/// @param t time
/// @param p path
/// @param dp derivative of path
/// @return lagranian evaluated at given time

template<typename T>
class kepler_lag{

  public :

  kepler_lag() {}

  T operator()(T t, 
  vector<T> p, vector<T> dp) const
  {
    if(p[0]<1e-2)
    return 0.5*(pow(dp[0], 2.0)+pow(p[0]*dp[1], 2.0))+100; // to avoid singularity  and penalize bad guess
    return 0.5*(pow(dp[0], 2.0)+pow(p[0]*dp[1], 2.0))+
    1/p[0];
  }

};

int main(void)
{
  PRECISION atol; // abs tol of action integral
  int num_eval; // number of points to eval
  int num_fourier; // number of sine and cosine function used for guess
  int max_iter; // number of iteration
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
  cout << " absolute tolerance of action: ";
  cin >> atol;
  cout << " number of sine and cosine function for approximation: ";
  cin >> num_fourier;
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
  #ifdef MONITOR
  cout << " file name to monitor optimization process: ";
  cin >> filename_monitor;
  #endif
  cout << " Now starts calculation" << endl;

  // initial condition
  libfourier::PI<PRECISION> pi;
  PRECISION zeta_max = zeta_min/(2*zeta_min-1);
  PRECISION a = (zeta_min+zeta_max)/2;
  PRECISION tmax = pi()*pow(a, 1.5);
  PRECISION period = 2*tmax;
  vector<PRECISION> p0 = {zeta_min, 0};
  vector<PRECISION> p1 = {zeta_max, pi()};
  mcm<PRECISION, kepler_lag<PRECISION>> kepler(PRECISION(0), 
  tmax, p0, p1, atol, num_fourier, period);

  kepler.set_init_guess();
  // variable to get optimization state
  int num_move; // number of actual moves
  PRECISION accept_ratio; // acceptance ratio

  #ifdef MONITOR
  tie(num_move, accept_ratio) = \
  kepler.optimize(max_iter, max_step, lambda, 
  filename_monitor, DIGITS);
  #else
  tie(num_move, accept_ratio) = \
  kepler.optimize(max_iter, max_step, lambda);
  #endif

  // Now fill time
  vector<PRECISION> t(num_eval+1, 0); // one more point need for end point
  PRECISION grid_space = tmax/PRECISION(num_eval);
  for(int i=1; i<num_eval; i++)
  t[i] = t[i-1]+grid_space;
  t[num_eval-1] = tmax;

  // store result
  vector<vector<PRECISION>> result(2, vector<PRECISION>(num_eval+1, 0));
  result = kepler.min_eval(t);
  min_action = kepler.get_min_action();

  // store coefficients
  vector<PRECISION> adder;
  vector<PRECISION> scaler;
  vector<vector<PRECISION>> coeff;

  tie(adder, scaler, coeff) = \
  kepler.get_min_coeff();

  int dim1 = adder.size();
  int dim2 = coeff[0].size();

  // move t by t0
  transform(t.begin(), t.end(), t.begin(),
  [t0](PRECISION &x){return x += t0;});

  cout << " Calcuation is finished" << endl;
  cout << "======================result==============================" << endl;
  cout.unsetf(ios::floatfield); // initialize floatfield
  cout.precision(DIGITS); // print significant digits
  cout << " Minimum action is " << min_action << endl;
  cout << " Number of Actual move is " << num_move << endl;
  cout << " Acceptance ratio: " << accept_ratio << endl;
  // store results to file
  fout.open(filename);
  fout << '#' << 't' << '\t' << "zeta" << '\t' << "theta" << endl;
  fout.unsetf(ios::floatfield); // initialize floatfield
  fout.precision(DIGITS); // print significant digits
  for(int i=0; i < num_eval+1; i++)
    fout << t[i] << '\t' << result[0][i] << '\t' << result[1][i] << endl;
  fout.close();
  fout.open(filename_coeff);
  fout.unsetf(ios::floatfield); // initialize floatfield
  fout.precision(DIGITS); // print significant digits
  fout << '#' << "adder" << '\t' << "scaler" << endl;
  for(int i=0; i<dim1; i++)
  fout << adder[i] << '\t' << scaler[i] << endl;
  fout << '#' << "zeta" << '\t' << "theta" << endl;
  for(int i=0; i<dim2; i++)
  fout << coeff[0][i] << '\t' << coeff[1][i] << endl;
  fout.close();
  cout << " Save result to " << filename << endl;
  cout << " Save coeffcients to " << filename_coeff << endl;
  cout << " Teriminates program, good bye :) " << endl;
  cout << "==========================================================" << endl;

  return 0;
}
      
      

  
						  
