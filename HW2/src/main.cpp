/*!
 * @file main.cpp
 * @brief main program for homework2 of Computer1 class in Yonsei University
 * Interactively reads inital condition, number of gird points to evaluate and
 * output file name then computes and saves solution.
 * @author pistack (Junho Lee)
 * @date 2021. 10. 28.
 */

#include <string>
#include <iostream>
#include <fstream>
#include "hw2.hpp"

using namespace std;

int main(void) {

  // init variables
  PRECISION zeta_min; // minimum value of zeta
  PRECISION t0; // initial time
  int n; // number of points to evaluate
  vector<PRECISION> t, y; // store results
  string filename; // file name to store results
  ofstream fout; // file output stream


  cout << "===============================================" << endl;
  cout << "         hw2: main program for homework2       " << endl;
  cout << " initial time: ";
  cin >> t0;
  cout << " minimum value of zeta: ";
  cin >> zeta_min;
  cout << " number of points to evaluate: ";
  cin >> n;
  cout << " file name to store result: ";
  cin >> filename;
  cout << " Now starts calculation" << endl;
  tie(t, y) = HW2<PRECISION>(zeta_min, t0, n);
  cout << " Calcuation is finished" << endl;

  // store result to file
  fout.open(filename);
  // first column t, second column zeta
  fout << '#' << 't' << '\t' << "zeta" << '\t' << endl;
  fout.unsetf(ios::floatfield); // initialize floatfield
  fout.precision(8); // print 8 significant digits
  for(int i=0; i < n+1; i++) // one more point needed for end point
      fout << t[i] << '\t' << y[i] << endl;
  fout.close();
  
  cout << " Save result to " << filename << endl;
  cout << " Teriminates program, good bye :) " << endl;
  cout << "===============================================" << endl;

  return 0;
}

  
     
			    
			  
    
