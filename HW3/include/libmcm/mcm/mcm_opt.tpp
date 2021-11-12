/*!
 * @file mcm_opt.tpp
 * @brief template for minimization of the action
 * using Monte Carlo Metropolis method
 * @author pistack (Junho Lee)
 * @date 2021. 11. 12.
 * @ingroup libmcm
 */

namespace libmcm {
template<typename T, typename Basis, typename Path, typename Lag>
std::tuple<std::size_t, T>
mcm<T, Basis, Path, Lag>::optimize(std::size_t num_iter, T step_size, T lambda)
{
  // number_of_accepted rand_walk
  std::size_t n_accept = 0;
  // acceptance ratio
  T prob = 0;

  // variable store temporal variable
  unsigned int dim_1 = init_path.size();
  T e; // estimated error
  T init_action; // initial action
  T accept_action; // accept_action
  T min_action; // minimum action
  std::vector<std::vector<T>> accept_guess = init_guess;
  std::vector<std::vector<T>> tmp_guess = init_guess;
  std::vector<std::vector<T>> min_guess = init_guess;
  std::vector<Path> tmp_path = init_path;

  mcm_action.update(init_path);
  init_action = mcm_action.eval(e);
  accept_action = init_action;
  min_action = init_action;

  for(std::size_t i=0; i<num_iter; ++i)
  {
    T r;
    T tmp_action;
    T delta_action;

    // sampling vaild path
    do
    {
      tmp_guess = rand_walk(accept_guess, step_size);
      for(unsigned int i=0; i<dim_1; ++i)
      {
        Basis tmp_basis(order, add_setup,
        tmp_guess[i]);
        tmp_path[i].update(tmp_basis);
      }
      mcm_action.update(tmp_path);
    }
    while(!mcm_action.is_vaild());

    // evaluate action
    tmp_action = mcm_action.eval(e);

    // accept rand_walk or not
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

template<typename T, typename Basis, typename Path, typename Lag>
std::tuple<std::size_t, T>
mcm<T, Basis, Path, Lag>::optimize(std::size_t num_iter, T step_size, T lambda,
std::string monitor)
{
  // number_of_accepted rand_walk
  std::size_t n_accept = 0;
  // acceptance ratio
  T prob = 0;

  // variable store temporal variable
  unsigned int dim_1 = init_path.size();
  T e; // estimated error
  T init_action; // initial action
  T accept_action; // accept_action
  T min_action; // minimum action
  std::vector<std::vector<T>> accept_guess = init_guess;
  std::vector<std::vector<T>> tmp_guess = init_guess;
  std::vector<std::vector<T>> min_guess = init_guess;
  std::vector<Path> tmp_path = init_path;

  mcm_action.update(init_path);
  init_action = mcm_action.eval(e);
  accept_action = init_action;
  min_action = init_action;

  // monitor optimization process
  std::ofstream fout(monitor, std::ios::out | std::ios::binary);

  for(int i=0; i<num_iter; ++i)
  {
    T r;
    T tmp_action;
    T delta_action;

    // sampling vaild path
    do
    {
      tmp_guess = rand_walk(accept_guess, step_size);
      for(int i=0; i<dim_1; ++i)
      {
        Basis tmp_basis(order, add_setup,
        tmp_guess[i]);
        tmp_path[i].update(tmp_basis);
      }
      mcm_action.update(tmp_path);
    }
    while(!mcm_action.is_vaild());

    // evaluate action
    tmp_action = mcm_action.eval(e);

    // accept rand_walk or not
    r = uniform_dist(gen);
    delta_action = tmp_action - accept_action;
    if(delta_action < 0 || r < std::exp(-lambda*delta_action))
    {
      n_accept++;
      accept_action = tmp_action;
      accept_guess = tmp_guess;
      // write accept_action to monitor
      struct W {
        T action, error;
      } w;
      w.action = accept_action;
      w.error = e;
      fout.write(reinterpret_cast<char *> (&w), sizeof(w));

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
}