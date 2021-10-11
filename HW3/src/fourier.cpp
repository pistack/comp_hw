/*!
 * @file fourier.cpp
 * @brief evaluate sum and derivative of sine and cosine function
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#include <cmath>
#include "basic.hpp"

using namespace std;

tuple<vector<double>, vector<double>>
sum_of_fourier(vector<double> &t, vector<double> c, int num_fourier)
{
  int max_sum = 2*num_fourier;
  int n = t.size();
  double half_period = t[n-1]-t[0];
  double omega = pi/half_period;
  double tmp = 0;
  vector<double> y(n, 0);
  vector<double> deriv(n, 0);

  for(int i = 0; i<n; i++)
    {
      for(int j = 0; j<max_sum; j += 2)
	{
	  tmp = (j/2+1)*omega;
	  if(c[j] != 0.0)
	    {
	      y[i] += c[j]*sin(tmp*t[i]);
	      deriv[i] += c[j]*tmp*cos(tmp*t[i]);
	    }
	  if(c[j+1] != 0.0)
	    {
	      y[i] += c[j+1]*cos(tmp*t[i]);
	      deriv[i] -= c[j+1]*tmp*sin(tmp*t[i]);
	    }
	}
    }

  return make_tuple(y, deriv);
}



