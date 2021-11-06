/*!
 * @file mcm_basic.tpp
 * @brief template which init, define initial guess
 * store minimum guess and evaluate guesses.
 * @author pistack (Junho Lee)
 * @date 2021. 11. 6.
 * @ingroup libmcm
 */

namespace libmcm {
template<typename T, typename Lag>
mcm<T, Lag> & mcm<T, Lag>::operator=(const mcm<T, Lag> &copy)
{
  t0 = copy.t0; t1 = copy.t1; p0 = copy.p0; p1 = copy.p1;
  num_fourier = copy.num_fourier; 
  fourier_period = copy.fourier_period;
  mcm_action = copy.mcm_action; init_guess = copy.init_guess;
  init_path = copy.init_path;
  min_guess = copy.min_guess; min_path = copy.min_path;
  return *this;
}

template<typename T, typename Lag>
void mcm<T, Lag>::set_init_guess(std::vector<std::vector<T>> init_c)
{
  int dim_1 = init_c.size();
  std::vector<libfourier::fourier_path<T>> \
  paths(dim_1, libfourier::fourier_path<T>());
  for(int i=0; i<dim_1; ++i)
  {
    libfourier::fourier<T> \
    tmp_fourier(num_fourier, fourier_period, init_c[i]);
    paths[i] = \
    libfourier::fourier_path<T>(t0, t1, p0[i], p1[i], tmp_fourier);
  }
  init_guess = init_c;
  init_path = paths;
}

template<typename T, typename Lag>
void mcm<T, Lag>::set_init_guess()
{
  int dim_1 = p0.size();
  int dim_2 = 2*num_fourier;
  std::vector<std::vector<T>> tmp_guess(dim_1, std::vector<T>(dim_2, 0));
  std::vector<libfourier::fourier_path<T>> paths(dim_1, 
  libfourier::fourier_path<T>());
  do
  {
    init_guess = move(tmp_guess, 1.0);
    for(int i=0; i<dim_1; ++i)
    {
      libfourier::fourier<T> \
      tmp_fourier(num_fourier, fourier_period, init_guess[i]);
      paths[i] = \
      libfourier::fourier_path<T>(t0, t1, p0[i], p1[i], tmp_fourier);
    }
    mcm_action.update(paths);
  } 
  while (! mcm_action.is_vaild());
  init_path = paths;
}

template<typename T, typename Lag>
T mcm<T, Lag>::get_init_action(T &e)
{
  T init_action;
  mcm_action.update(init_path);
  init_action = mcm_action.eval(e);
  return init_action;
}

template<typename T, typename Lag>
std::tuple<std::vector<T>, std::vector<T>, std::vector<std::vector<T>>>
mcm<T, Lag>::get_init_coeff()
{
  int dim_1 = init_path.size();
  std::vector<T> init_adder(dim_1, 0);
  std::vector<T> init_scaler(dim_1, 0);
  for(int i=0; i<dim_1; ++i)
  {
    init_adder[i] = init_path[i].get_adder();
    init_scaler[i] = init_path[i].get_scaler();
  }
  return std::make_tuple(init_adder, init_scaler, init_guess);
}

template<typename T, typename Lag>
std::vector<T> mcm<T, Lag>::init_eval(T t)
{
  int dim_1 = init_path.size();
  std::vector<T> result(dim_1, 0);
  for(int i=0; i<dim_1; ++i)
  result[i] = init_path[i].eval(t);
  return result;
}

template<typename T, typename Lag>
std::vector<std::vector<T>> mcm<T, Lag>::init_eval(std::vector<T> t)
{
  int dim_1 = init_path.size();
  int dim_2 = t.size();
  std::vector<std::vector<T>> result(dim_1, std::vector<T>(dim_2, 0));
  for(int i=0; i<dim_1; ++i)
  result[i] = init_path[i].eval(t);
  return result;
}

template<typename T, typename Lag>
T mcm<T, Lag>::get_min_action(T &e)
{
  T min_action;
  mcm_action.update(min_path);
  min_action = mcm_action.eval(e);
  return min_action;
}

template<typename T, typename Lag>
std::tuple<std::vector<T>, std::vector<T>, std::vector<std::vector<T>>>
mcm<T, Lag>::get_min_coeff()
{
  int dim_1 = min_path.size();
  std::vector<T> min_adder(dim_1, 0);
  std::vector<T> min_scaler(dim_1, 0);
  for(int i=0; i<dim_1; ++i)
  {
    min_adder[i] = min_path[i].get_adder();
    min_scaler[i] = min_path[i].get_scaler();
  }
  return std::make_tuple(min_adder, min_scaler, min_guess);
}

template<typename T, typename Lag>
std::vector<T> mcm<T, Lag>::min_eval(T t)
{
  int dim_1 = min_path.size();
  std::vector<T> result(dim_1, 0);
  for(int i=0; i<dim_1; ++i)
  result[i] = min_path[i].eval(t);
  return result;
}

template<typename T, typename Lag>
std::vector<std::vector<T>> mcm<T, Lag>::min_eval(std::vector<T> t)
{
  int dim_1 = min_path.size();
  int dim_2 = t.size();
  std::vector<std::vector<T>> result(dim_1, std::vector<T>(dim_2, 0));
  for(int i=0; i<dim_1; ++i)
  result[i] = min_path[i].eval(t);
  return result;
}
}