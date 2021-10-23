/*!
 * @file main.cpp
 * @brief main program for homework3 of Computer1 class in Yonsei University
 * Interactively reads inital condition, number of sine function used for guess,
 * number of gird points to evaluate, number of interation, step size and
 * output file name then computes and saves solution.
 * @author pistack (Junho Lee)
 * @date 2021. 10. 24.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include "mcm.hpp"

const double pi = 3.141592653589793; ///<define pi

using namespace std;

int main(void)
{
  double atol, rtol; // abs and rel tol of action integral
  int num_eval; // number of points to eval
  int num_fourier; // number of sine and cosine function used for guess
  int max_iter; // number of iteration
  double zeta_min; // minimum value of zeta
  double t0; // initial time
  double max_step; // step size
  double lambda; // parameter to accept move
  double min_action; // minimum action value
  string filename; // file name to store results
  string filename_coeff; // file name to store coeffcients
  string filename_monitor; // file name to monitor optimization process
  ofstream fout; // file output stream

  cout << "==========================================================" << endl;
  cout << "               hw3: main program for homework3            " << endl;
  cout << " initial time: ";
  cin >> t0;
  cout << " minimum value of zeta: ";
  cin >> zeta_min;
  cout << " absolute tolerance of action: ";
  cin >> atol;
  cout << " relative tolerance of action: ";
  cin >> rtol;
  cout << " number of sine and cosine function for approximation: ";
  cin >> num_fourier;
  cout << " number of points to evaluate path: ";
  cin >> num_eval;
  cout << " size of step: ";
  cin >> max_step;
  cout << " value for paramter to accept move: ";
  cin >> lambda;
  cout << " Number of iteration: ";
  cin >> max_iter;
  cout << " file name to store result: ";
  cin >> filename;
  cout << " file name to store coefficient: ";
  cin >> filename_coeff;
  cout << " file name to monitor optimization process: ";
  cin >> filename_monitor;
  cout << " Now starts calculation" << endl;

  // initial condition
  double zeta_max = zeta_min/(2*zeta_min-1);
  double a = 0.5*(zeta_min+zeta_max);
  double tmax = pi*pow(a, 1.5);
  double period = 2*tmax;
  vector<double> p0 = {zeta_min, 0.0};
  vector<double> p1 = {zeta_max, pi};
  mcm kepler(0.0, tmax, p0, p1, atol, rtol, num_fourier,
  period, [](double t, vector<double> x, vector<double> dx)
  { if(abs(x[0])>1e-8)
    return 0.5*(pow(dx[0],2.0)+pow(x[0]*dx[1], 2.0))+1/abs(x[0]);
    else
    return 0.0;});

  kepler.set_init_guess();
  
  
  if(num_fourier > 1)
  {
    vector<vector<double>> guess_coeff(2, vector<double>(2*num_fourier, 0));
    guess_coeff[0][0] = -0.71876708;
    guess_coeff[0][1] = -0.53458926;
    guess_coeff[1][0] = -0.027881983;
    guess_coeff[1][1] = 0.28787967;
    kepler.set_init_guess(guess_coeff);
  }

  if(num_fourier > 2)
  {
    vector<vector<double>> guess_coeff(2, vector<double>(2*num_fourier, 0));
    guess_coeff[0][0] = -0.57920731;
    guess_coeff[0][1] = -0.89324102;
    guess_coeff[0][2] = -0.077075367;
    guess_coeff[0][3] = -0.36377004;
    guess_coeff[1][0] = -0.12627861;
    guess_coeff[1][1] = 0.99027267;
    guess_coeff[1][2] = -0.22496388;
    guess_coeff[1][3] = 0.016038242;
    kepler.set_init_guess(guess_coeff);
  } 
  
  // variable to get optimization state
  int num_move; // number of actual moves
  double accept_ratio; // acceptance ratio

  tie(num_move, accept_ratio) = \
  kepler.optimize(max_iter, max_step, lambda, filename_monitor);

  // Now fill time
  vector<double> t(num_eval, 0);
  double grid_space = tmax/double(num_eval-1);
  for(int i=1; i<num_eval; i++)
  t[i] = t[i-1]+grid_space;
  t[num_eval-1] = tmax;

  // store result
  vector<vector<double>> result(2, vector<double>(num_eval, 0));
  result = kepler.min_eval(t);
  min_action = kepler.get_min_action();

  // store coefficients
  vector<double> adder;
  vector<double> scaler;
  vector<vector<double>> coeff;

  tie(adder, scaler, coeff) = \
  kepler.get_min_coeff();

  int dim1 = adder.size();
  int dim2 = coeff[0].size();

  // move t by t0
  transform(t.begin(), t.end(), t.begin(),
  [t0](double &x){return x += t0;});

  cout << " Calcuation is finished" << endl;
  cout << "======================result==============================" << endl;
  cout.unsetf(ios::floatfield); // initialize floatfield
  cout.precision(8); // print 8 significant digits
  cout << " Minimum action is " << min_action << endl;
  cout << " Number of Actual move is " << num_move << endl;
  cout << " Acceptance ratio: " << accept_ratio << endl;
  // store results to file
  fout.open(filename);
  fout << '#' << 't' << '\t' << "zeta" << '\t' << "theta" << endl;
  fout.unsetf(ios::floatfield); // initialize floatfield
  fout.precision(8); // print 8 significant digits
  for(int i=0; i < num_eval; i++)
    fout << t[i] << '\t' << result[0][i] << '\t' << result[1][i] << endl;
  fout.close();
  fout.open(filename_coeff);
  fout.unsetf(ios::floatfield); // initialize floatfield
  fout.precision(8); // print 8 significant digits
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
      
      

  
						  
