/*!
 * @file main.cpp
 * @brief main program for homework3 of Computer1 class in Yonsei University
 * Interactively reads inital condition, number of sine function used for guess,
 * number of gird points to evaluate, number of interation, step size and
 * output file name then computes and saves solution.
 * @author pistack (Junho Lee)
 * @date 2021. 10. 10.
 */

#include <string>
#include <iostream>
#include <fstream>
#include "hw3.hpp"

using namespace std;

int main(void)
{
  // init variables
  random_device rd;
  mt19937 gen(rd()); // set random number generator
  uniform_real_distribution<double> dist(-1, 1);   // set distribution
  int num_action; // number of points to approximate kepler action integral
  int num_eval; // number of points to eval
  int num_fourier; // number of sine and cosine function used for guess
  int num_iter; // number of iteration
  int num_move; // number of actual moves
  double zeta_min; // minimum value of zeta
  double t0; // initial time
  double step; // step size
  double lambda; // parameter for adapt step size
  double min_action; // minimum action value
  vector<double> t(num_eval,0); // time
  vector<double> zeta(num_eval, 0); // zeta
  vector<double> theta(num_eval, 0); // theta
  string filename; // file name to store results
  ofstream fout; // file output stream

  cout << "==========================================================" << endl;
  cout << "               hw3: main program for homework3            " << endl;
  cout << " initial time: ";
  cin >> t0;
  cout << " minimum value of zeta: ";
  cin >> zeta_min;
  cout << " number of points used for action: ";
  cin >> num_action;
  cout << " number of sine and cosine function for approximation: ";
  cin >> num_fourier;
  cout << " number of points to evaluate path: ";
  cin >> num_eval;
  cout << " size of step: ";
  cin >> step;
  cout << " value for paramter to adapt step size: ";
  cin >> lambda;
  cout << " number of iteration: ";
  cin >> num_iter;
  cout << " file name to store result: ";
  cin >> filename;
  cout << " Now starts calculation" << endl;
  tie(num_move, min_action, t, zeta, theta) =				\
    HW3(t0, zeta_min, 
    num_action, num_fourier, 
    num_eval, num_iter, step, lambda, gen, dist);
  cout << " Calcuation is finished" << endl;
  cout << "======================result==============================" << endl;
  cout << " Minimum action is " << min_action << endl;
  cout << " Number of Actual move is " << num_move << endl;

  // store results to file
  fout.open(filename);
  fout << '#' << 't' << '\t' << "zeta" << '\t' << "theta" << endl;
  fout.unsetf(ios::floatfield); // initialize floatfield
  fout.precision(8); // print 8 significant digits
  for(int i=0; i < num_eval; i++)
    fout << t[i] << '\t' << zeta[i] << '\t' << theta[i] << endl;
  fout.close();
  
  cout << " Save result to " << filename << endl;
  cout << " Teriminates program, good bye :) " << endl;
  cout << "==========================================================" << endl;

  return 0;
}
      
      

  
						  
