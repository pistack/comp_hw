/*!
 * @file kepler.cpp
 * @brief evaluates action for kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#include <cmath>
#include "kepler.hpp"

using namespace std;


double eval_action(vector<double> t,
fourier_path zeta, fourier_path theta)
{
  int n = t.size();
  double action = 0;
  vector<double> t_mid(n-1, 0);
  vector<double> zeta_t(n, 0);
  vector<double> deriv_zeta_t(n, 0);
  vector<double> zeta_t_mid(n-1, 0);
  vector<double> deriv_zeta_t_mid(n-1, 0);
  vector<double> deriv_theta_t(n, 0);
  vector<double> deriv_theta_t_mid(n-1, 0);
 

  // fill mid points
  for(int i=1; i<n; i++)
    t_mid[i-1] = 0.5*(t[i-1]+t[i]);

  // evaluate path
  zeta_t = zeta.eval(t);
  zeta_t_mid = zeta.eval(t_mid);
  deriv_zeta_t = zeta.deriv(t);
  deriv_zeta_t_mid = zeta.deriv(t_mid);
  deriv_theta_t = theta.deriv(t);
  deriv_theta_t_mid = theta.deriv(t_mid);

  
  // Evaulate action by simple's 1/3 rule
  
  for(int i=1; i<n; i++)
    {
      double dt = t[i]-t[i-1];
      double tmp1;
      double tmp2;
      double tmp3;

      // initial point
      tmp1 = 0.5*(pow(deriv_zeta_t[i-1], 2.0) + \
		  pow((zeta_t[i-1]*deriv_theta_t[i-1]), 2.0)) + \
      1/abs(zeta_t[i-1]);

      // mid point
      tmp2 = 0.5*(pow(deriv_zeta_t_mid[i-1], 2.0) + \
		  pow((zeta_t_mid[i-1]*deriv_theta_t_mid[i-1]), 2.0)) + \
      1/abs(zeta_t_mid[i-1]);

      // end point
      tmp3 = 0.5*(pow(deriv_zeta_t[i], 2.0) +		    \
		  pow((zeta_t[i]*deriv_theta_t[i]), 2.0)) + \
      1/abs(zeta_t[i]);
      
      action += dt/6.0*(tmp1+4.0*tmp2+tmp3);
    }
  return action;
  
}
