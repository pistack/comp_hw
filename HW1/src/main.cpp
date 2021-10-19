/*!
 * @file main.cpp
 * @brief main program for homework1 of Computer1 class in Yonsei University
 * Interactively reads inital condition, number of gird points to evaluate and
 * output file name then computes and saves solution.
 * @author pistack (Junho Lee)
 * @date 2021. 10. 15.
 */

#include <string>
#include <iostream>
#include <fstream>
#include "hw1.hpp"

using namespace std;

int main(void) {

  // init variables
  double t0; // initial time
  double t1; // finial time
  int n; // number of points to evaluate
  double y0; // intial condition zeta(0)
  double y0p; // initial condition zeta'(0)
  double theta0; // initial condition theta(0)
  vector<double> t, y, theta; // store results
  string filename; // file name to store results
  ofstream fout; // file output stream
  


  cout << "===============================================" << endl;
  cout << "         hw1: main program for homework1       " << endl;
  cout << " initial time: ";
  cin >> t0;
  cout << "  finial time: ";
  cin >> t1;
  cout << " number of points to evaluate: ";
  cin >> n;
  cout << " initial value of zeta: ";
  cin >> y0;
  cout << " initial value of the derivative of zeta: ";
  cin >> y0p;
  cout << " initial value of theta: ";
  cin >> theta0;
  cout << " file name to store result: ";
  cin >> filename;
  cout << " Now starts calculation" << endl;
  tie(t, y, theta) = HW1(t0, t1, n, y0, y0p, theta0);
  cout << " Calcuation is finished" << endl;

  // store result to file
  fout.open(filename);
  // first column t, second column zeta
  fout << '#' << 't' << '\t' << "zeta" << '\t' \
  << "theta" << endl;
  fout.unsetf(ios::floatfield); // initialize floatfield
  fout.precision(8); // print 8 significant digits
  for(int i=0; i < n+1; i++) // one more point needed for end point
      fout << t[i] << '\t' << y[i] << '\t' << theta[i] << endl;
  fout.close();
  
  cout << " Save result to " << filename << endl;
  cout << " Teriminates program, good bye :) " << endl;
  cout << "===============================================" << endl;

  return 0;
}
