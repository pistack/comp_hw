/*!
 * @file main.cpp
 * @brief main program for homework3 of Computer1 class in Yonsei University
 * Interactively reads inital condition, number of sine function used for guess,
 * number of gird points to evaluate, number of interation, step size and
 * output file name then computes and saves solution.
 * @author pistack (Junho Lee)
 * @date 2021. 10. 10.
 */

#include <algorithm>
#include <cmath>
#include <cerrno>
#include <string>
#include <iostream>
#include <fstream>
#include "hw3.hpp"

const double pi = 3.141592653589793; ///<define pi

using namespace std;

int main(void)
{
  double atol, rtol; // abs and rel tol of action integral
  int num_eval; // number of points to eval
  int num_fourier; // number of sine and cosine function used for guess
  int max_iter; // number of iteration
  double conv_prob; // converge criteria 
  // abs and rel change of minimum action
  double zeta_min; // minimum value of zeta
  double t0; // initial time
  double max_step; // step size
  double lambda; // parameter for adapt step size
  double min_action; // minimum action value
  string filename; // file name to store results
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
  cout << " value for paramter to adapt step size: ";
  cin >> lambda;
  cout << " Maximum number of iteration: ";
  cin >> max_iter;
  cout << " Converge criteria: probability: ";
  cin >> conv_prob;
  cout << " file name to store result: ";
  cin >> filename;
  cout << " Now starts calculation" << endl;

  // initial condition
  double zeta_max = zeta_min/(2*zeta_min-1);
  double a = 0.5*(zeta_min+zeta_max);
  double tmax = pi*pow(a, 1.5);
  double period = 2*tmax;
  vector<double> p0 = {zeta_min, 0.0};
  vector<double> p1 = {zeta_max, pi};
  HW3 kepler(0.0, tmax, p0, p1, atol, rtol, num_fourier,
  period, [](double t, vector<double> x, vector<double> dx)
  {return 0.5*(pow(dx[0],2.0)+pow(x[0]*dx[1], 2.0))+1/abs(x[0]);});

  // initial guess
  // 2 of dimension of path
  // 2*num_fourier for the number of terms
  vector<vector<double>> c(2, vector<double>(2*num_fourier, 0));
  c[0][0] = 0.5; c[0][1] = -0.5; c[1][0] = 0.5; c[1][1] = -0.5;

  kepler.set_init_guess(c);
  // variable to get optimization state
  int num_move; // number of actual moves
  int num_conv; // number of iteration to converge
  double opt_prob; // probability that guess is minimum

  tie(num_move, num_conv, opt_prob) = \
  kepler.optimize(max_iter, max_step, 
  lambda, conv_prob);

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

  // move t by t0
  transform(t.begin(), t.end(), t.begin(),
  [t0](double &x){return x += t0;});

  cout << " Calcuation is finished" << endl;
  cout << "======================result==============================" << endl;
  if(errno == ERANGE)
  cout << " Iteration does not meet converge ceriteria!" << endl;
  else
  cout << " Convergence criteria meets at " << num_conv << " iterations!" << endl;
  cout.unsetf(ios::floatfield); // initialize floatfield
  cout.precision(8); // print 8 significant digits
  cout << " Minimum action is " << min_action << endl;
  cout << " Number of Actual move is " << num_move << endl;
  cout << " Probability that guess is minimum: " << opt_prob << endl;
  // store results to file
  fout.open(filename);
  fout << '#' << 't' << '\t' << "zeta" << '\t' << "theta" << endl;
  fout.unsetf(ios::floatfield); // initialize floatfield
  fout.precision(8); // print 8 significant digits
  for(int i=0; i < num_eval; i++)
    fout << t[i] << '\t' << result[0][i] << '\t' << result[1][i] << endl;
  fout.close();
  cout << " Save result to " << filename << endl;
  cout << " Teriminates program, good bye :) " << endl;
  cout << "==========================================================" << endl;

  return 0;
}
      
      

  
						  
