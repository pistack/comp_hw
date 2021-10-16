/*!
 * @file hw3.cpp
 * @brief code for homework3 of Computer1 class in Yonsei University
 *        Minimize the action by Markov Chain Monte Carlo Method 
 *        to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#include <algorithm>
#include <cmath>
#include "fourier_path.hpp"
#include "action.hpp"
#include "hw3.hpp"

using namespace std;

tuple<int, double, vector<double>, vector<double>, vector<double>>
HW3(double t0, double zeta_min, double atol, double rtol, 
int num_fourier, int num_eval,
int num_iter, double step, double lambda,
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
  double period = 2*tmax;

  // variable used for adaptation of step size
  double delta_action;
  double adapt_step = step;

  // variable used to store minimal action
  double min_action;

  // guess minimial
  vector<double> min_c_zeta(size, 0);
  vector<double> min_c_theta(size, 0);

  // store path
  vector<double> path_zeta_t(num_eval, 0);
  vector<double> path_theta_t(num_eval, 0);

  // initial guess
  min_c_zeta[0] = 0.5;
  min_c_zeta[1] = -0.5;
  min_c_theta[0] = 0.5;
  min_c_theta[1] = -0.5;

  // init fourier function
  fourier fourier_zeta(num_fourier, period, min_c_zeta);
  fourier fourier_theta(num_fourier, period, min_c_theta);

  // init path
  vector<fourier_path> path;
  path.push_back(fourier_path(0.0, tmax, zeta_min, zeta_max, fourier_zeta));
  path.push_back(fourier_path(0.0, tmax, 0.0, pi, fourier_theta));

  // init action
  action kepler_action(atol, rtol, path,
  [](double t, vector<double> p, vector<double> dp)
  {return 0.5*(pow(dp[0], 2.0)+pow(p[0]*dp[1], 2.0))+1/abs(p[0]);});

  min_action = kepler_action.eval();

  for(int i = 0; i < num_iter; i++)
    {
      // init temporal variable
      double tmp_action;

      // guess temporal
      vector<double> tmp_c_zeta(size, 0);
      vector<double> tmp_c_theta(size, 0);

      do
      {
        tmp_c_zeta = move_step(min_c_zeta, step, gen, dist);
        tmp_c_theta = move_step(min_c_theta, step, gen, dist);
        // update fourier coefficients
        fourier_zeta.update(tmp_c_zeta);
        fourier_theta.update(tmp_c_theta);
        path[0].update(fourier_zeta);
        path[1].update(fourier_theta);
        kepler_action.update(path);
      }
      while(!kepler_action.is_vaild());

      tmp_action = kepler_action.eval();

      delta_action = tmp_action - min_action;
      
      // adapt step size
      if(adapt_step > exp(-lambda*delta_action*adapt_step))
      adapt_step = exp(-lambda*delta_action*adapt_step);

      // update action and guess
      if(delta_action < 0)
      {
        min_action = tmp_action;
        min_c_zeta = tmp_c_zeta;
        min_c_theta = tmp_c_theta;
        n_accept++; // movement is accepted when action decreases
      }
    }

    // update fourier coeffcients to minimum one
    fourier_zeta.update(min_c_zeta);
    fourier_theta.update(min_c_theta);
    
    // update fourier function to minimum one
    path[0].update(fourier_zeta);
    path[1].update(fourier_theta);
    
    // minimum action for debugging
    // kepler_action.update(path_zeta, path_theta);
    // min_action = kepler_action.eval();

    // initialize time
    vector<double> t(num_eval, 0);
    // fill time
    for(int i=1; i<num_eval; i++)
    t[i] = tmax*double(i)/double(num_eval-1);
    path_zeta_t = path[0].eval(t);
    path_theta_t = path[1].eval(t);

    transform(t.begin(), t.end(), t.begin(),
    [t0](double &x){return x += t0;});

  return make_tuple(n_accept, min_action, t,
  path_zeta_t, path_theta_t);
}
