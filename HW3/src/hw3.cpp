/*!
 * @file hw3.cpp
 * @brief code for homework3 of Computer1 class in Yonsei University
 *        Minimize the action by Markov Chain Monte Carlo Method 
 *        to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 10.
 */

#include <iostream>
#include "hw3.hpp"

using namespace std;

tuple<int, double, vector<double>, vector<double>, vector<double>>
HW3(double zeta_min, double t0, int n, int num_fourier, int num_iter,
    double step, double lambda,
    mt19937 &gen, uniform_real_distribution<double> &dist)
{
  int n_accept = 0;
  int size = 2*num_fourier;
  double zeta_max = zeta_min/(2*zeta_min-1);
  double a = 0.5*(zeta_min+zeta_max);
  double tmax = pi*pow(a,1.5);
  double delta_action;
  double adapt_step = step;
  double scale_zeta;
  double scale_theta;
  double add_zeta;
  double add_theta;
  double min_action;
  vector<double> t(n, 0);
  vector<double> min_c_zeta(size, 0);
  vector<double> min_c_theta(size, 0);
  vector<double> tmp_c_zeta(size, 0);
  vector<double> tmp_c_theta(size, 0);
  vector<double> min_zeta(n, 0);
  vector<double> min_theta(n, 0);
  vector<double> min_deriv_zeta(n, 0);
  vector<double> min_deriv_theta(n, 0);

  for(int i = 0; i<n; i++)
    {
      t[i] = double(i)/double(n-1)*tmax;
    }

  min_c_zeta[0] = 0.5;
  min_c_zeta[1] = -0.5;
  min_c_theta[0] = 0.5;
  min_c_theta[1] = -0.5;

  tie(min_zeta, min_deriv_zeta) = \
    sum_of_fourier(t, min_c_zeta, num_fourier);
  tie(min_theta, min_deriv_theta) = \
    sum_of_fourier(t, min_c_theta, num_fourier);
  
  scale_zeta = (zeta_max-zeta_min)/(min_zeta[n-1]-min_zeta[0]);
  add_zeta = zeta_min - scale_zeta*min_zeta[0];
  scale_theta = pi/(min_theta[n-1]-min_theta[0]);
  add_theta = -scale_theta*min_theta[0];
  
  min_zeta = scale_and_add_vector(min_zeta, scale_zeta, add_zeta);
  min_theta = scale_and_add_vector(min_theta, scale_theta, add_theta);
  
  min_deriv_zeta = scale_and_add_vector(min_deriv_zeta, scale_zeta, 0);
  min_deriv_theta = scale_and_add_vector(min_deriv_theta, scale_theta, 0);
  
  min_action = eval_action(t, min_zeta, min_deriv_zeta,
			   min_theta, min_deriv_theta);

  for(int i = 0; i < num_iter; i++)
    {
      double scale_zeta, scale_theta;
      double tmp_action;
      vector<double> tmp_zeta(n, 0);
      vector<double> tmp_theta(n, 0);
      vector<double> tmp_deriv_zeta(n, 0);
      vector<double> tmp_deriv_theta(n, 0);
      do
	{
	  tmp_c_zeta = move_step(min_c_zeta, step, gen, dist);
          tmp_c_theta = move_step(min_c_theta, step, gen, dist);
          tie(tmp_zeta, tmp_deriv_zeta) = \
	    sum_of_fourier(t, tmp_c_zeta, num_fourier);
	  tie(tmp_theta, tmp_deriv_theta) = \
	    sum_of_fourier(t, tmp_c_theta, num_fourier);
	}
      while(abs(tmp_zeta[0]-tmp_zeta[n-1]) < 1e-8 ||
	    abs(tmp_theta[0]-tmp_theta[n-1]) < 1e-8);

      scale_zeta = (zeta_max-zeta_min)/(tmp_zeta[n-1]-tmp_zeta[0]);
      add_zeta = zeta_min - scale_zeta*tmp_zeta[0];
      scale_theta = pi/(tmp_theta[n-1]-tmp_theta[0]);
      add_theta = -scale_theta*tmp_theta[0];
      
      tmp_zeta = scale_and_add_vector(tmp_zeta, scale_zeta, add_zeta);
      tmp_deriv_zeta = scale_and_add_vector(tmp_deriv_zeta, scale_zeta, 0);
      
      tmp_theta = scale_and_add_vector(tmp_theta, scale_theta, add_theta);
      tmp_deriv_theta = scale_and_add_vector(tmp_deriv_theta, scale_theta, 0);
      
      tmp_action = eval_action(t, tmp_zeta, tmp_deriv_zeta,
			       tmp_theta, tmp_deriv_theta);

      delta_action = tmp_action - min_action;
      if(adapt_step > exp(-lambda*delta_action*adapt_step))
	adapt_step = exp(-lambda*delta_action*adapt_step);

      if(tmp_action < min_action)
	{
	  min_action = tmp_action;
	  min_c_zeta = tmp_c_zeta;
	  min_c_theta = tmp_c_theta;
	  min_zeta = tmp_zeta;
	  min_theta = tmp_theta;
	  n_accept++;
	}
    }

  return make_tuple(n_accept, min_action,
		    scale_and_add_vector(t,1, t0),
		    min_zeta,
		    min_theta);
}
