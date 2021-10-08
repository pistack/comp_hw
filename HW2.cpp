/*
HW2
C++ code for homework2 in computer1 class in Yonsei University
Using numerical integration to solve Kepler problem
*/

#include <cmath>
#include <tuple>
#include <vector>
#include <iostream>

using namespace std;

tuple<vector<double>, vector<double>> HW2(double a,
					  double t0,
					  double zeta_min, double zeta_max,
					  int n)
{

  //initialize variable
  double grid_space;
  vector<double> zeta(n+1, 0);
  vector<double> t(n+1, 0);
  vector<double> x1(n+1, 0);
  vector<double> x2(n+1, 0);

  // use uniform n pts between zeta_min and zeta_max
  // additional one point needed for end point.

  grid_space = (zeta_max-zeta_min)/n;
  
  for(int i=0; i<n+1; i++)
    {
      zeta[i] = i*grid_space + zeta_min;
    }

  // change of variable
  // and seperate integral
  for(int i=0; i<n+1; i++)
    {
      x1[i] = sqrt(zeta[i]-zeta_min);
      x2[i] = sqrt(zeta_max-zeta[i]);
    }

  // initial condtion zeta(t0) = zeta_min
  t[0] = t0;

  // numerical integration using Riemann-Stieltiges integral
  // Quite fun, integrand_1(x1) = x2 and
  // integrand_2(x2) = x1
  for(int i=1; i<n+1; i++)
    {
      t[i] = t[i-1] + 8/sqrt(a)*(x1[i]-x1[i-1])*x2[i-1] - \
	10/sqrt(a)*(x2[i]-x2[i-1])*x1[i-1];
    }

  return make_tuple(t, zeta);
}

int main(void) {
  double zeta_min = 0.9;
  double zeta_max = 9./8.;
  int n = 20;
  double t0 = 0;
  double a = 0.987654;
  vector<double> t, y;
  
  tie(t, y) = HW2(a, t0, zeta_min, zeta_max, n);

    for (int i = 0; i < n+1; i++)
      {
        cout << t[i] << '\t' << y[i] << endl;
      }

    return 0;
}

  
     
			    
			  
    
