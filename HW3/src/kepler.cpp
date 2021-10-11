/*!
 * @file kepler.cpp
 * @brief evaluates action for kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#include <cmath>
#include "basic.hpp"

using namespace std;


double eval_action(vector<double> &t,
		   double zeta_init, double zeta_finial,
		   double theta_init, double theta_finial,
		   vector<double> c_zeta,
		   vector<double> c_theta)
{
  int n = t.size();
  vector<double> t_mid(n+1, 0);
  vector<double> zeta(n, 0);
  vector<double> deriv_zeta(n, 0);
  vector<double> zeta_mid(n+1, 0);
  vector<double> deriv_zeta_mid(n+1, 0);
  vector<double> theta(n, 0);
  vector<double> deriv_theta(n, 0);
  vector<double> theta_mid(n+1, 0);
  vector<double> deriv_theta_mid(n+1, 0);
  double action = 0;

  // fill mid points
  t_mid[0] = t[0];
  t_mid[n] = t[n-1];
  for(int i=1; i<n; i++)
    t_mid[i] = 0.5*(t[i-1]+t[i]);

  // evaluate path

  // zeta
  tie(zeta, deriv_zeta) = eval_path(t, zeta_init, zeta_finial, c_zeta);
  tie(zeta_mid, deriv_zeta_mid) = \
    eval_path(t_mid, zeta_init, zeta_finial, c_zeta);

  // theta
  tie(theta, deriv_theta) = eval_path(t, theta_init, theta_finial, c_theta);
  tie(theta_mid, deriv_theta_mid) = \
    eval_path(t_mid, theta_init, theta_finial, c_theta);

  
  // Evaulate action by simple's 1/3 rule
  
  for(int i=1; i<n; i++)
    {
      double dt = t[i]-t[i-1];
      double tmp1;
      double tmp2;
      double tmp3;

      // initial point
      tmp1 = 0.5*(pow(deriv_zeta[i-1], 2.0) + \
		  pow((zeta[i-1]*deriv_theta[i-1]), 2.0)) + \
	1/abs(zeta[i-1]);

      // mid point
      tmp2 = 0.5*(pow(deriv_zeta_mid[i], 2.0) + \
		  pow((zeta_mid[i]*deriv_theta_mid[i]), 2.0)) + \
	1/abs(zeta_mid[i]);

      // end point
      tmp3 = 0.5*(pow(deriv_zeta[i], 2.0) +		    \
		  pow((zeta[i]*deriv_theta[i]), 2.0)) + \
	1/abs(zeta[i]);
      
      action += dt/6.0*(tmp1+4.0*tmp2+tmp3);
    }
  return action;
  
}
