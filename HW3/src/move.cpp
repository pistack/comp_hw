/*!
 * @file move.cpp
 * @brief randomly move initial guess, at most step size.
 * @author pistack (Junho Lee)
 * @date 2021. 10. 18.
 * @ingroup hw3
 */

#include <cmath>
#include <numeric>
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
      moved_guess[i][j] = guess[i][j] + max_step*normal_dist(gen);
      // bound guess to [-1, 1]
      if(moved_guess[i][j]>1)
      moved_guess[i][j] = 1;
      if(moved_guess[i][j]<-1)
      moved_guess[i][j] = -1;
    }
  }
  return moved_guess;
}
