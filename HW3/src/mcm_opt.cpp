/*!
 * @file mcm_opt.cpp
 * @brief code to minimize the action
 * using Monte Carlo Metropolis method
 * @author pistack (Junho Lee)
 * @date 2021. 10. 24.
 * @ingroup libmcm
 */

#include <cmath>
#include <fstream>
#include "mcm.hpp"

using namespace std;

std::tuple<int, double>
mcm::optimize(int num_iter, double step_size, double lambda)
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

  for(int i=0; i<num_iter; i++)
  {
    double r;
    double tmp_action;
    double delta_action;

    // sampling vaild path
    do
    {
      tmp_guess = move(accept_guess, step_size);
      for(int i=0; i<dim_1; i++)
      {
        fourier tmp_fourier(num_fourier, fourier_period,
        tmp_guess[i]);
        tmp_path[i].update(tmp_fourier);
      }
      mcm_action.update(tmp_path);
    }
    while(!mcm_action.is_vaild());

    // evaluate action
    tmp_action = mcm_action.eval();

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
    prob = double(n_accept)/double(num_iter);
  }
  return make_tuple(n_accept, prob);
}

std::tuple<int, double>
mcm::optimize(int num_iter, double step_size, double lambda,
string monitor)
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

  // monitor optimization process
  ofstream fout;
  fout.open(monitor);
  fout.unsetf(ios::floatfield);
  fout.precision(8);


  for(int i=0; i<num_iter; i++)
  {
    double r;
    double tmp_action;
    double delta_action;

    // sampling vaild path
    do
    {
      tmp_guess = move(accept_guess, step_size);
      for(int i=0; i<dim_1; i++)
      {
        fourier tmp_fourier(num_fourier, fourier_period,
        tmp_guess[i]);
        tmp_path[i].update(tmp_fourier);
      }
      mcm_action.update(tmp_path);
    }
    while(!mcm_action.is_vaild());

    // evaluate action
    tmp_action = mcm_action.eval();

    // accept move or not
    r = uniform_dist(gen);
    delta_action = tmp_action - accept_action;
    if(delta_action < 0 || r < exp(-lambda*delta_action))
    {
      n_accept++;
      accept_action = tmp_action;
      accept_guess = tmp_guess;
      // write accept_action to monitor
      fout << n_accept << '\t' << accept_action << endl;
      if(accept_action < min_action)
      {
        min_action = accept_action;
        min_guess = accept_guess;
        min_path = tmp_path; // in that case tmp_path == accept_path
      }
    }
    prob = double(n_accept)/double(num_iter);
  }
  // file should be closed
  fout.close();
  return make_tuple(n_accept, prob);
}