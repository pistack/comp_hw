/*!
 * @file hw2.cpp
 * @brief code for homework2 of Computer1 class in Yonsei University
 * Use numerical integration to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 10.
 */

#include "hw2.hpp"

using namespace std;

tuple<vector<double>, vector<double>>
HW2(double zeta_min, double t0, int n)
{

  //initialize variable
  double grid_space; // grid spacing
  double a; 
  double zeta_max; // maximum value of zeta
  double tmp;
  double c1, c2;
  vector<double> zeta(n+1, 0);
  vector<double> t(n+1, 0);
  vector<double> x1(n+1, 0);
  vector<double> x2(n+1, 0);

  // Note a = zeta_min**(-2)-2*zeta_min**(-1)
  //        = zeta_max**(-2)-2*zeta_max**(-1)

  tmp = 1/zeta_min;
  a = tmp*(tmp-2);
  tmp = 1.0+sqrt(1.0+a);
  zeta_max = -tmp/a;

  // calculates c1 and c2
  tmp = zeta_max - zeta_min;
  tmp = 1/tmp;
  c1 = 2*zeta_min*tmp;
  c2 = 2*zeta_max*tmp;
  tmp = sqrt(-a);
  tmp = 1/tmp;
  c1 *= tmp;
  c2 *= tmp;

  // use uniform n pts between zeta_min and zeta_max
  // additional one point needed for end point.

  grid_space = (zeta_max-zeta_min)/double(n);
  
  for(int i=0; i<n+1; i++)
    {
      zeta[i] = double(i)*grid_space + zeta_min;
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
      t[i] = t[i-1] + c1*(x1[i]-x1[i-1])*x2[i-1] - \
	c2*(x2[i]-x2[i-1])*x1[i-1];
    }

  return make_tuple(t, zeta);
}
