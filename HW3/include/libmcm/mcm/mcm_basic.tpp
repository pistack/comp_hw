/*!
 * @file mcm_basic.tpp
 * @brief template which init, define initial guess
 * store minimum guess and evaluate guesses.
 * @author pistack (Junho Lee)
 * @date 2021. 11. 12.
 * @ingroup libmcm
 */

namespace libmcm {
template<typename T, typename Basis, typename Path, typename Lag>
mcm<T, Basis, Path, Lag> & mcm<T, Basis, Path, Lag>::operator=(const mcm<T, Basis, Path, Lag> &copy)
{
  t0 = copy.t0; t1 = copy.t1; p0 = copy.p0; p1 = copy.p1;
  order = copy.order; 
  add_setup = copy.add_setup;
  mcm_action = copy.mcm_action;
  init_path = copy.init_path;
  min_path = copy.min_path;
  return *this;
}

template<typename T, typename Basis, typename Path, typename Lag>
void mcm<T, Basis, Path, Lag>::set_init_guess(std::vector<std::vector<T>> init_c)
{
  unsigned int dim_1 = init_c.size();
  std::vector<Path> \
  paths(dim_1, Path());
  for(unsigned int i=0; i<dim_1; ++i)
  {
    Basis \
    tmp_basis(order, add_setup, init_c[i]);
    paths[i] = \
    Path(t0, t1, p0[i], p1[i], tmp_basis);
  }
  init_guess = init_c;
  init_path = paths;
}

template<typename T, typename Basis, typename Path, typename Lag>
void mcm<T, Basis, Path, Lag>::set_init_guess()
{
  unsigned int dim_1 = p0.size(); // dimension of path
  unsigned int dim_2 = Basis().terms(order); // number of terms
  // tmp guess
  std::vector<std::vector<T>> tmp_guess(dim_1, std::vector<T>(dim_2, 0));
  std::vector<Path> paths(dim_1, Path());
  do
  {
    init_guess = rand_walk(tmp_guess, 1.0); // sample initial guess randomly
    for(unsigned int i=0; i<dim_1; ++i)
    {
      Basis \
      tmp_basis(order, add_setup, init_guess[i]);
      paths[i] = \
      Path(t0, t1, p0[i], p1[i], tmp_basis);
    }
    mcm_action.update(paths);
  } 
  while (! mcm_action.is_vaild());
  init_path = paths;
}

template<typename T, typename Basis, typename Path, typename Lag>
T mcm<T, Basis, Path, Lag>::get_init_action(T &e)
{
  T init_action;
  mcm_action.update(init_path);
  init_action = mcm_action.eval(e);
  return init_action;
}

template<typename T, typename Basis, typename Path, typename Lag>
std::vector<T> mcm<T, Basis, Path, Lag>::init_eval(T t)
{
  unsigned int dim_1 = init_path.size();
  std::vector<T> result(dim_1, 0);
  for(unsigned int i=0; i<dim_1; ++i)
  result[i] = init_path[i].eval(t);
  return result;
}

template<typename T, typename Basis, typename Path, typename Lag>
std::vector<std::vector<T>> mcm<T, Basis, Path, Lag>::init_eval(std::vector<T> t)
{
  unsigned int dim_1 = init_path.size();
  std::size_t dim_2 = t.size();
  std::vector<std::vector<T>> result(dim_1, std::vector<T>(dim_2, 0));
  for(unsigned int i=0; i<dim_1; ++i)
  result[i] = init_path[i].eval(t);
  return result;
}

template<typename T, typename Basis, typename Path, typename Lag>
T mcm<T, Basis, Path, Lag>::get_min_action(T &e)
{
  T min_action;
  mcm_action.update(min_path);
  min_action = mcm_action.eval(e);
  return min_action;
}

template<typename T, typename Basis, typename Path, typename Lag>
std::vector<T> mcm<T, Basis, Path, Lag>::min_eval(T t)
{
  unsigned int dim_1 = min_path.size();
  std::vector<T> result(dim_1, 0);
  for(unsigned int i=0; i<dim_1; ++i)
  result[i] = min_path[i].eval(t);
  return result;
}

template<typename T, typename Basis, typename Path, typename Lag>
std::vector<std::vector<T>> mcm<T, Basis, Path, Lag>::min_eval(std::vector<T> t)
{
  unsigned int dim_1 = min_path.size();
  std::size_t dim_2 = t.size();
  std::vector<std::vector<T>> result(dim_1, std::vector<T>(dim_2, 0));
  for(unsigned int i=0; i<dim_1; ++i)
  result[i] = min_path[i].eval(t);
  return result;
}
}