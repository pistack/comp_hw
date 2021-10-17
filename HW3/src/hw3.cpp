/*!
 * @file hw3.cpp
 * @brief code for homework3 of Computer1 class in Yonsei University
 *        Minimize the action by Markov Chain Monte Carlo Method 
 *        to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 * @ingroup hw3
 */

#include <algorithm>
#include <cmath>
#include <cerrno>
#include <numeric>
#include "fourier.hpp"
#include "fourier_path.hpp"
#include "action.hpp"
#include "hw3.hpp"

using namespace std;

double HW3::dist(std::vector<std::vector<double>> x,
std::vector<std::vector<double>> y)
{
  int dim_1 = x.size();
  int dim_2 = x[0].size();
  vector<double> z(dim_2, 0);
  vector<double> norms(dim_1, 0);

  for(int i=0; i<dim_1; i++)
  {
    for(int j=0; j<dim_2; j++)
    z[j] = x[i][j] - y[i][j];
    norms[i] = sqrt(inner_product(z.begin(), z.end(), 
    z.begin(), 0.0));
  }
  return sqrt(inner_product(norms.begin(), norms.end(), 
  norms.begin(), 0.0));
}

std::vector<std::vector<double>>
HW3::move(std::vector<std::vector<double>> guess, double max_step)
{
  int dim_1 = guess.size();
  int dim_2 = guess[0].size();
  vector<vector<double>> moved_guess(dim_1,
  vector<double>(dim_2, 0));

  for(int i=0; i<dim_1; i++)
  {
    for(int j=0; j<dim_2; j++)
    {
      moved_guess[i][j] = guess[i][j] + max_step*uniform_dist(gen);
      // bound guess to [-1, 1]
      if(moved_guess[i][j]>1)
      moved_guess[i][j] = 1;
      if(moved_guess[i][j]<-1)
      moved_guess[i][j] = -1;
    }
  }
  return moved_guess;
}

HW3 & HW3::operator=(const HW3 &copy)
{
  t0 = copy.t0; t1 = copy.t1; p0 = copy.p0; p1 = copy.p1;
  num_fourier = copy.num_fourier; hw3_period = copy.hw3_period;
  hw3_action = copy.hw3_action; init_guess = copy.init_guess;
  init_path = copy.init_path; init_action = copy.init_action;
  min_guess = copy.min_guess; min_path = copy.min_path;
  min_action = copy.min_action;
  return *this;
}

void HW3::set_init_guess(std::vector<vector<double>> init_c)
{
  int dim_1 = init_c.size();
  vector<fourier_path> paths(dim_1, fourier_path());
  for(int i=0; i<dim_1; i++)
  {
    fourier tmp_fourier(num_fourier, hw3_period, init_c[i]);
    paths[i] = fourier_path(t0, t1, p0[i], p1[i], tmp_fourier);
  }
  init_guess = init_c;
  init_path = paths;
  hw3_action.update(init_path);
  init_action = hw3_action.eval();
  min_guess = init_guess;
  min_path = init_path;
  min_action = init_action;
}

double HW3::get_init_action()
{
  return init_action;
}

tuple<vector<double>, vector<double>, vector<vector<double>>>
HW3::get_init_coeff()
{
  int dim_1 = init_path.size();
  vector<double> init_adder(dim_1, 0);
  vector<double> init_scaler(dim_1, 0);
  for(int i=0; i<dim_1; i++)
  {
    init_adder[i] = init_path[i].get_adder();
    init_scaler[i] = init_path[i].get_scaler();
  }
  return make_tuple(init_adder, init_scaler, init_guess);
}

vector<double> HW3::init_eval(double t)
{
  int dim_1 = init_path.size();
  vector<double> result(dim_1, 0);
  for(int i=0; i<dim_1; i++)
  result[i] = init_path[i].eval(t);
  return result;
}

vector<vector<double>> HW3::init_eval(vector<double> t)
{
  int dim_1 = init_path.size();
  int dim_2 = t.size();
  vector<vector<double>> result(dim_1, vector<double>(dim_2, 0));
  for(int i=0; i<dim_1; i++)
  result[i] = init_path[i].eval(t);
  return result;
}

double HW3::get_min_action()
{
  return min_action;
}

tuple<vector<double>, vector<double>, vector<vector<double>>>
HW3::get_min_coeff()
{
  int dim_1 = min_path.size();
  vector<double> min_adder(dim_1, 0);
  vector<double> min_scaler(dim_1, 0);
  for(int i=0; i<dim_1; i++)
  {
    min_adder[i] = min_path[i].get_adder();
    min_scaler[i] = min_path[i].get_scaler();
  }
  return make_tuple(min_adder, min_scaler, min_guess);
}

vector<double> HW3::min_eval(double t)
{
  int dim_1 = min_path.size();
  vector<double> result(dim_1, 0);
  for(int i=0; i<dim_1; i++)
  result[i] = min_path[i].eval(t);
  return result;
}

vector<vector<double>> HW3::min_eval(vector<double> t)
{
  int dim_1 = min_path.size();
  int dim_2 = t.size();
  vector<vector<double>> result(dim_1, vector<double>(dim_2, 0));
  for(int i=0; i<dim_1; i++)
  result[i] = min_path[i].eval(t);
  return result;
}

std::tuple<int, int, double>
HW3::optimize(int max_iter, double max_step, double lambda,
double conv_prob)
{
  // number_of_accepted move
  int n_accept = 0;
  // probability that guess is minimum
  double prob = 0;

  // check_converge
  int n_conv = 0;
  bool is_converged = false;

  // variable used to adapt step size

  double adapt_step = max_step;

  // variable store temporal variable
  int dim_1 = init_path.size();
  vector<vector<double>> tmp_guess = init_guess;
  vector<fourier_path> tmp_path = init_path;

  for(int i=0; i<max_iter; i++)
  {
    double tmp_action;
    double delta_action;

    n_conv=i+1;
    // sampling vaild path
    do
    {
      tmp_guess = move(min_guess, adapt_step);
      for(int i=0; i<dim_1; i++)
      {
        fourier tmp_fourier(num_fourier, hw3_period,
        tmp_guess[i]);
        tmp_path[i].update(tmp_fourier);
      }
      hw3_action.update(tmp_path);
    }
    while(!hw3_action.is_vaild());

    // evaluate action
    tmp_action = hw3_action.eval();

    // adapt step
    delta_action = tmp_action - min_action;
    if(adapt_step > exp(-lambda*delta_action*adapt_step))
    adapt_step = exp(-lambda*delta_action*adapt_step);

    // update min_action, min_guess, min_path
    // if action decreases
    // (i.e. move is accepted)
    if(delta_action < 0)
    {
      n_accept++;
      // now update
      min_action = tmp_action;
      min_guess = tmp_guess;
      min_path = tmp_path;
      // reset adapt step
      adapt_step = max_step;
    }
    // check convergence
    prob = double(n_conv-n_accept-1)/double(n_conv);
    if(prob>conv_prob)
      {
        is_converged = true;
        break;
      }
  }
  if(!is_converged)
  errno=ERANGE;   
  return make_tuple(n_accept, n_conv, prob);
}