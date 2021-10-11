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
  int n; // number of points to evaluate
  int num_fourier; // number of sine and cosine function used for guess
  int num_iter; // number of iteration
  int num_move; // number of actual moves
  double zeta_min; // minimum value of zeta
  double t0; // initial time
  double step; // step size
  double lambda; // parameter for adapt step size
  double min_action; // minimum action value
  vector<double> t(n,0); // time
  vector<double> zeta(n,0);  // guessed zeta
  vector<double> theta(n,0); // guessed theta
  string filename_zeta; // file name to store results zeta
  string filename_theta; // file name to store results theta
  ofstream fout; // file output stream

  cout << "===============================================" << endl;
  cout << "         hw3: main program for homework3       " << endl;
  cout << " initial time: ";
  cin >> t0;
  cout << " minimum value of zeta: ";
  cin >> zeta_min;
  cout << " number of points to evaluate: ";
  cin >> n;
  cout << " number of sine and cosine function used for guess: ";
  cin >> num_fourier;
  cout << " size of step: ";
  cin >> step;
  cout << " value for paramter to adapt step size: ";
  cin >> lambda;
  cout << " number of iteration: ";
  cin >> num_iter;
  cout << " file name to store result (zeta): ";
  cin >> filename_zeta;
  cout << " file name to store result (theta): ";
  cin >> filename_theta;
  cout << " Now starts calculation" << endl;
  tie(num_move, min_action, t, zeta, theta) =				\
    HW3(zeta_min, t0, n, num_fourier, num_iter, step, lambda, gen, dist);
  cout << " Calcuation is finished" << endl;
  cout << "======================result===================" << endl;
  cout << " Minimum action is " << min_action << endl;
  cout << " Number of Actual move is " << num_move << endl;
  fout.open(filename_zeta);
  for(int i=0; i < n; i++) 
    {
      fout << t[i] << '\t' << zeta[i] << endl;
    }
  fout.close();
  cout << " Save result (zeta) to " << filename_zeta << endl;
  fout.open(filename_theta);
  for(int i=0; i < n; i++) 
    {
      fout << t[i] << '\t' << theta[i] << endl;
    }
  fout.close();
  cout << " Save result (theta) to " << filename_theta << endl;
  cout << " Teriminates program, good bye :) " << endl;
  cout << "===============================================" << endl;

  return 0;
}
      
      

  
						  
