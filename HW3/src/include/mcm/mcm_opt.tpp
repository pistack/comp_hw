/*!
 * @file mcm_opt.tpp
 * @brief template for minimization of the action
 * using Monte Carlo Metropolis method
 * @author pistack (Junho Lee)
 * @date 2021. 10. 28.
 * @ingroup libmcm
 */

template<typename T>
std::tuple<int, T>
mcm<T>::optimize(int num_iter, T step_size, T lambda)
{
  // number_of_accepted move
  int n_accept = 0;
  // acceptance ratio
  T prob = 0;

  // variable store temporal variable
  int dim_1 = init_path.size();
  T accept_action = init_action;
  std::vector<std::vector<T>> accept_guess = init_guess;
  std::vector<std::vector<T>> tmp_guess = init_guess;
  std::vector<fourier_path<T>> tmp_path = init_path;

  for(int i=0; i<num_iter; i++)
  {
    T r;
    T tmp_action;
    T delta_action;

    // sampling vaild path
    do
    {
      tmp_guess = move(accept_guess, step_size);
      for(int i=0; i<dim_1; i++)
      {
        fourier<T> tmp_fourier(num_fourier, fourier_period,
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
    if(delta_action < 0 || r < std::exp(-lambda*delta_action))
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
    prob = T(n_accept)/T(num_iter);
  }
  return std::make_tuple(n_accept, prob);
}

template<typename T>
std::tuple<int, T>
mcm<T>::optimize(int num_iter, T step_size, T lambda,
std::string monitor)
{
  // number_of_accepted move
  int n_accept = 0;
  // acceptance ratio
  T prob = 0;

  // variable store temporal variable
  int dim_1 = init_path.size();
  T accept_action = init_action;
  std::vector<std::vector<T>> accept_guess = init_guess;
  std::vector<std::vector<T>> tmp_guess = init_guess;
  std::vector<fourier_path<T>> tmp_path = init_path;

  // monitor optimization process
  std::ofstream fout;
  fout.open(monitor);
  fout.unsetf(std::ios::floatfield);
  fout.precision(8);


  for(int i=0; i<num_iter; i++)
  {
    T r;
    T tmp_action;
    T delta_action;

    // sampling vaild path
    do
    {
      tmp_guess = move(accept_guess, step_size);
      for(int i=0; i<dim_1; i++)
      {
        fourier<T> tmp_fourier(num_fourier, fourier_period,
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
    if(delta_action < 0 || r < std::exp(-lambda*delta_action))
    {
      n_accept++;
      accept_action = tmp_action;
      accept_guess = tmp_guess;
      // write accept_action to monitor
      fout << n_accept << '\t' << accept_action << std::endl;
      if(accept_action < min_action)
      {
        min_action = accept_action;
        min_guess = accept_guess;
        min_path = tmp_path; // in that case tmp_path == accept_path
      }
    }
    prob = T(n_accept)/T(num_iter);
  }
  // file should be closed
  fout.close();
  return std::make_tuple(n_accept, prob);
}