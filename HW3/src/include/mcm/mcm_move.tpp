/*!
 * @file mcm_move.tpp
 * @brief sample guess via ramdom walk.
 * @author pistack (Junho Lee)
 * @date 2021. 10. 28.
 * @ingroup libmcm
 */

template<typename T>
std::vector<std::vector<T>>
mcm<T>::move(std::vector<std::vector<T>> guess, T step_size)
{
  int dim_1 = guess.size();
  int dim_2 = guess[0].size();
  std::vector<std::vector<T>> moved_guess(dim_1,
  std::vector<T>(dim_2, 0));

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
