/*!
 * @file hw3.cpp
 * @brief code for homework3 of Computer1 class in Yonsei University
 * Minimize the action by Monte Carlo Metropolis
 * method described in 
 * [Entropy 2020, 22(9), 916](https://doi.org/10.3390/e22090916)
 * @author pistack (Junho Lee)
 * @date 2021. 10. 19.
 * @ingroup hw3
 */

#include <cmath>
#include "fourier.hpp"
#include "fourier_path.hpp"
#include "action.hpp"
#include "hw3.hpp"

using namespace std;

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

void HW3::set_init_guess()
{
  int dim_1 = p0.size();
  int dim_2 = 2*num_fourier;
  vector<vector<double>> tmp_guess(dim_1, vector<double>(dim_2, 0));
  vector<fourier_path> paths(dim_1, fourier_path());
  do
  {
    init_guess = move(tmp_guess, 1.0);
    for(int i=0; i<dim_1; i++)
    {
      fourier tmp_fourier(num_fourier, hw3_period, init_guess[i]);
      paths[i] = fourier_path(t0, t1, p0[i], p1[i], tmp_fourier);
    }
    hw3_action.update(paths);
  } 
  while (! hw3_action.is_vaild());
  init_path = paths;
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

std::tuple<int, double>
HW3::optimize(int max_iter, double max_step, double lambda)
{
  // number_of_accepted move
  int n_accept = 0;
  // acceptance ratio
  double prob = 0;

  // variable store temporal variable
  int dim_1 = init_path.size();
  double accept_action = init_action;
  vector<vector<double>> accept_guess = init_guess;
  vector<vector<double>> tmp_guess = init_guess;
  vector<fourier_path> tmp_path = init_path;

  for(int i=0; i<max_iter; i++)
  {
    double r;
    double tmp_action;
    double delta_action;

    // sampling vaild path
    do
    {
      tmp_guess = move(accept_guess, max_step);
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

    // accept move or not
    r = uniform_dist(gen);
    delta_action = tmp_action - accept_action;
    if(delta_action < 0 || r < exp(-lambda*delta_action))
    {
      n_accept++;
      accept_action = tmp_action;
      accept_guess = tmp_guess;
      if(accept_action < min_action)
      {
        min_action = accept_action;
        min_guess = accept_guess;
        min_path = tmp_path; // in that case tmp_path == accept_path
      }
    }
    prob = double(n_accept)/double(max_iter);
  }
  return make_tuple(n_accept, prob);
}