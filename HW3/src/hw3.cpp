/*!
 * @file hw3.cpp
 * @brief code for homework3 of Computer1 class in Yonsei University
 *        Minimize the action by Markov Chain Monte Carlo Method 
 *        to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#include <cmath>
#include "basic.hpp"
#include "kepler.hpp"

using namespace std;

tuple<int, double, vector<double>, vector<double>, vector<double>>
HW3(double zeta_min, double t0, int n, int num_fourier, int num_iter,
    double step, double lambda,
    mt19937 &gen, uniform_real_distribution<double> &dist)
{
  // number of accepted move
  int n_accept = 0;

  // number of terms (one fourier = one sine + one cosine => 2 terms)
  int size = 2*num_fourier;

  // inital condition
  double zeta_max = zeta_min/(2*zeta_min-1);
  double a = 0.5*(zeta_min+zeta_max);
  double tmax = pi*pow(a,1.5);

  // variable used for adaptation of step size
  double delta_action;
  double adapt_step = step;

  // variable used to store minimal action
  double min_action;

  // time
  vector<double> t(n, 0);
  
  // subset of time to check our guess is vailed
  vector<double> tst(2, 0);

  // guess minimial
  vector<double> min_c_zeta(size, 0);
  vector<double> min_c_theta(size, 0);

  // guess temporal
  vector<double> tmp_c_zeta(size, 0);
  vector<double> tmp_c_theta(size, 0);

  // minimal path
  vector<double> min_zeta(n, 0);
  vector<double> min_theta(n, 0);
  vector<double> min_deriv_zeta(n, 0);
  vector<double> min_deriv_theta(n, 0);

  // fill time
  for(int i = 0; i<n; i++)
    t[i] = double(i)/double(n-1)*tmax;
  
  tst[0] = t[0];
  tst[1] = t[n-1];

  // initial guess
  min_c_zeta[0] = 0.5;
  min_c_zeta[1] = -0.5;
  min_c_theta[0] = 0.5;
  min_c_theta[1] = -0.5;

  // initial action
  min_action = eval_action(t, zeta_min, zeta_max, 0.0, pi,
			   min_c_zeta, min_c_theta);

  for(int i = 0; i < num_iter; i++)
    {
      double tmp_action;
      vector<double> tmp_zeta(2, 0);
      vector<double> tmp_theta(2, 0);
      vector<double> tmp_deriv_zeta(2, 0);
      vector<double> tmp_deriv_theta(2, 0);
      do
	{
	  tmp_c_zeta = move_step(min_c_zeta, step, gen, dist);
          tmp_c_theta = move_step(min_c_theta, step, gen, dist);
          tie(tmp_zeta, tmp_deriv_zeta) = \
	    sum_of_fourier(tst, tmp_c_zeta, num_fourier);
	  tie(tmp_theta, tmp_deriv_theta) = \
	    sum_of_fourier(tst, tmp_c_theta, num_fourier);
	}
      while(abs(tmp_zeta[0]-tmp_zeta[1]) < 1e-8 ||
	    abs(tmp_theta[0]-tmp_theta[1]) < 1e-8);
      
      tmp_action = eval_action(t, zeta_min, zeta_max, 0.0, pi,
			       tmp_c_zeta, tmp_c_theta);

      delta_action = tmp_action - min_action;
      if(adapt_step > exp(-lambda*delta_action*adapt_step))
	adapt_step = exp(-lambda*delta_action*adapt_step);

      if(tmp_action < min_action)
	{
	  min_action = tmp_action;
	  min_c_zeta = tmp_c_zeta;
	  min_c_theta = tmp_c_theta;
	  n_accept++;
	}
    }

  // evaluate minimal path
  tie(min_zeta, min_deriv_zeta) = eval_path(t, zeta_min, zeta_max,
					    min_c_zeta);
  tie(min_theta, min_deriv_theta) = eval_path(t, 0.0, pi,
					      min_c_theta);

  return make_tuple(n_accept, min_action,
		    scale_and_add_vector(t,1, t0),
		    min_zeta,
		    min_theta);
}
