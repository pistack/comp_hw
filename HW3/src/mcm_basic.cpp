/*!
 * @file mcm_basic.cpp
 * @brief code to init, define initial guess
 * store minimum guess and evaluate guesses.
 * @author pistack (Junho Lee)
 * @date 2021. 10. 24.
 * @ingroup libmcm
 */

#include "mcm.hpp"

using namespace std;

mcm & mcm::operator=(const mcm &copy)
{
  t0 = copy.t0; t1 = copy.t1; p0 = copy.p0; p1 = copy.p1;
  num_fourier = copy.num_fourier; 
  fourier_period = copy.fourier_period;
  mcm_action = copy.mcm_action; init_guess = copy.init_guess;
  init_path = copy.init_path; init_action = copy.init_action;
  min_guess = copy.min_guess; min_path = copy.min_path;
  min_action = copy.min_action;
  return *this;
}

void mcm::set_init_guess(std::vector<vector<double>> init_c)
{
  int dim_1 = init_c.size();
  vector<fourier_path> paths(dim_1, fourier_path());
  for(int i=0; i<dim_1; i++)
  {
    fourier tmp_fourier(num_fourier, fourier_period, init_c[i]);
    paths[i] = fourier_path(t0, t1, p0[i], p1[i], tmp_fourier);
  }
  init_guess = init_c;
  init_path = paths;
  mcm_action.update(init_path);
  init_action = mcm_action.eval();
  min_guess = init_guess;
  min_path = init_path;
  min_action = init_action;
}

void mcm::set_init_guess()
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
      fourier tmp_fourier(num_fourier, fourier_period, init_guess[i]);
      paths[i] = fourier_path(t0, t1, p0[i], p1[i], tmp_fourier);
    }
    mcm_action.update(paths);
  } 
  while (! mcm_action.is_vaild());
  init_path = paths;
  init_action = mcm_action.eval();
  min_guess = init_guess;
  min_path = init_path;
  min_action = init_action;
}

double mcm::get_init_action()
{
  return init_action;
}

tuple<vector<double>, vector<double>, vector<vector<double>>>
mcm::get_init_coeff()
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

vector<double> mcm::init_eval(double t)
{
  int dim_1 = init_path.size();
  vector<double> result(dim_1, 0);
  for(int i=0; i<dim_1; i++)
  result[i] = init_path[i].eval(t);
  return result;
}

vector<vector<double>> mcm::init_eval(vector<double> t)
{
  int dim_1 = init_path.size();
  int dim_2 = t.size();
  vector<vector<double>> result(dim_1, vector<double>(dim_2, 0));
  for(int i=0; i<dim_1; i++)
  result[i] = init_path[i].eval(t);
  return result;
}

double mcm::get_min_action()
{
  return min_action;
}

tuple<vector<double>, vector<double>, vector<vector<double>>>
mcm::get_min_coeff()
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

vector<double> mcm::min_eval(double t)
{
  int dim_1 = min_path.size();
  vector<double> result(dim_1, 0);
  for(int i=0; i<dim_1; i++)
  result[i] = min_path[i].eval(t);
  return result;
}

vector<vector<double>> mcm::min_eval(vector<double> t)
{
  int dim_1 = min_path.size();
  int dim_2 = t.size();
  vector<vector<double>> result(dim_1, vector<double>(dim_2, 0));
  for(int i=0; i<dim_1; i++)
  result[i] = min_path[i].eval(t);
  return result;
}