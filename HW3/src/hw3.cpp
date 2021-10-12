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

  // period of fourier function
  double period;

  // variable used for adaptation of step size
  double delta_action;
  double adapt_step = step;

  // variable used to store minimal action
  double min_action;

  // time
  vector<double> t(n, 0);
  
  // subset of time to check our guess is vailed
  vector<double> tends(2, 0) = {0.0, tmax};

  // guess minimial
  vector<double> min_c_zeta(size, 0);
  vector<double> min_c_theta(size, 0);

  // guess temporal
  vector<double> tmp_c_zeta(size, 0);
  vector<double> tmp_c_theta(size, 0);

  // minimal path
  fourier_path min_zeta = fourier_path();
  fourier_path min_theta = fourier_path();

  // fill time
  for(int i = 0; i<n; i++)
    t[i] = double(i)/double(n-1)*tmax;

  // set period of fourier function
  period = 2*tmax;

  // initial guess
  min_c_zeta[0] = 0.5;
  min_c_zeta[1] = -0.5;
  min_c_theta[0] = 0.5;
  min_c_theta[1] = -0.5;

  // init path
  min_zeta.init(0.0, tmax, zeta_min, zeta_max,
  period, min_c_zeta);
  min_theta.init(0.0, tmax, 0.0, pi,
  period, min_c_theta);

  // initial action
  min_action = eval_action(t, zeta_min, min_zeta, min_theta);

  for(int i = 0; i < num_iter; i++)
    {
      // init temporal variable
      double tmp_action;
      fourier_path tmp_zeta = fourier_path();
      fourier_path tmp_theta = fourier_path();

      do
      {
        // init temporal variable to check guess is vaild
        fourier tmp_zeta_check = fourier();
        fourier tmp_theta_check = fourier();
        vector<double> tst_zeta(2, 0);
        vector<double> tst_theta(2, 0);

        tmp_c_zeta = move_step(min_c_zeta, step, gen, dist);
        tmp_c_theta = move_step(min_c_theta, step, gen, dist);
        tmp_zeta_check.init(num_fourier, period,
        tmp_c_zeta);
        tmp_theta_check.init(num_fourier, period,
        tmp_c_theta);
        tst_zeta = tmp_theta_check.eval(tends);
        tst_theta = tmp_theta_check.eval(tends);
      }
      while(abs(tst_zeta[0]-tst_zeta[1]) < 1e-8 ||
      abs(tst_theta[0]-tst_theta[1]) < 1e-8);
      
      tmp_zeta.init(0.0, tmax, zeta_min, zeta_max, period,
      tmp_c_zeta);
      tmp_theta.init(0.0, tmax, zeta_min, zeta_max, period,
      tmp_c_theta);
      tmp_action = eval_action(t, tmp_zeta, tmp_theta);
      
      delta_action = tmp_action - min_action;
      
      // adapt step size
      if(adapt_step > exp(-lambda*delta_action*adapt_step))
      adapt_step = exp(-lambda*delta_action*adapt_step);

      // update action and guess
      if(tmp_action < min_action)
      {
        min_action = tmp_action;
        min_c_zeta = tmp_c_zeta;
        min_c_theta = tmp_c_theta;
        n_accept++; // move is accepted when action decreases
      }
    }

  // minimal path
  min_zeta.init(0.0, tmax, zeta_min, zeta_max,
  period, min_c_zeta);
  min_theta.init(0.0, tmax, 0.0, pi,
  period, min_c_theta);

  // minimal action
  min_action = eval_action(t, zeta_min, min_zeta, min_theta);

  return make_tuple(n_accept, min_action,
		    scale_and_add_vector(t, 1.0, t0),
		    min_zeta.eval(t),
		    min_theta.eval(t));
}
