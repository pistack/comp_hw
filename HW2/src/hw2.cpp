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
  vector<double> u(n+1, 0);
  vector<double> v(n+1, 0);
  vector<double> u_mid(n, 0);
  vector<double> v_mid(n, 0);

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
  zeta[i] = double(i)*grid_space + zeta_min;

  // change of variable
  // and seperate integral
  for(int i=1; i<n; i++)
    {
      u[i] = sqrt(zeta[i]-zeta_min);
      v[i] = sqrt(zeta_max-zeta[i]);
      u_mid[i-1] = sqrt(0.5*(zeta[i]+zeta[i-1])-zeta_min);
      v_mid[i-1] = sqrt(zeta_max-0.5*(zeta[i]+zeta[i-1]));
    }
  // end point
  tmp = sqrt(zeta_max-zeta_min);
  u[n] = tmp;
  v[0] = tmp;
  u_mid[n-1] = sqrt(0.5*(zeta_max+zeta[n-1])-zeta_min);
  v_mid[n-1] = sqrt(0.5*(zeta_max-zeta[n-1]));
  

  // initial condtion zeta(t0) = zeta_min
  t[0] = t0;

  // numerical integration using 
  // non equi-spacing Simpson's rule
  // Quite fun, integrand_1(x1) = x2 and
  // integrand_2(x2) = x1
  for(int i=1; i<n+1; i++)
    {
      double h0 = u_mid[i-1]-u[i-1];
      double h1 = u[i] - u_mid[i-1];
      double d = u[i] - u[i-1];

      t[i] = t[i-1];

      // u part
      tmp = d/6.0*((2.0-h1/h0)*v[i-1]+
      (d/h0)*(d/h1)*v_mid[i-1]+
      (2.0-h0/h1)*v[i]);
      t[i] += c1*tmp;

      // v part
      h0 = v_mid[i-1]-v[i-1];
      h1 = v[i] - v_mid[i-1];
      d = v[i]-v[i-1];
      tmp = d/6.0*((2.0-h1/h0)*u[i-1]+
      (d/h0)*(d/h1)*u_mid[i-1]+
      (2.0-h0/h1)*u[i]);
      t[i] -= c2*tmp;
    }

  return make_tuple(t, zeta);
}
