/*!
 * @file mcm_move.cpp
 * @brief sample guess via ramdom walk.
 * @author pistack (Junho Lee)
 * @date 2021. 10. 24.
 * @ingroup libmcm
 */

#include "mcm.hpp"

using namespace std;

std::vector<std::vector<double>>
mcm::move(std::vector<std::vector<double>> guess, double step_size)
{
  int dim_1 = guess.size();
  int dim_2 = guess[0].size();
  vector<vector<double>> moved_guess(dim_1,
  vector<double>(dim_2, 0));

  for(int i=0; i<dim_1; i++)
  {
    for(int j=0; j<dim_2; j++)
    {
      moved_guess[i][j] = guess[i][j] + step_size*normal_dist(gen);
      // bound guess to [-1, 1]
      if(moved_guess[i][j]>1)
      moved_guess[i][j] = 1;
      if(moved_guess[i][j]<-1)
      moved_guess[i][j] = -1;
    }
  }
  return moved_guess;
}
