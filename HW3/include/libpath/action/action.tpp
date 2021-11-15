/*!
 * @file action.tpp
 * @ingroup libpath
 * @brief evaluates action
 * @author pistack (Junho Lee)
 * @date 2021. 11. 15.
 */

namespace libpath {

/// check vaildity 
template<typename T, typename Path, typename Lag>
void action<T, Path, Lag>::check_vaild()
{
  unsigned int n = path_action.size();
  vaildity = false;

  for(unsigned int i=0; i<n; ++i)
  {
    if(! path_action[i].is_vaild())
    return;
  }
  vaildity = true;
  return;
}

template<typename T, typename Path, typename Lag>
T action<T, Path, Lag>::eval_lagrangian(T t)
{
  unsigned int n = path_action.size();
  std::vector<T> p(n, 0); // path
  std::vector<T> dp(n, 0); // derivative of path

  for(unsigned int i=0; i<n; ++i)
  {
    p[i] = path_action[i].eval(t);
    dp[i] = path_action[i].deriv(t);
  }
  return lag(t, p, dp);
}

template<typename T, typename Path, typename Lag>
std::vector<T> action<T, Path, Lag>::eval_lagrangian(std::vector<T> t)
{
  unsigned int n = path_action.size();
  unsigned int m = t.size();
  std::vector<T> p(n, 0); // path
  std::vector<T> dp(n, 0); // derivative of path
  std::vector<T> lag_val(m, 0); // result

  for(unsigned int j=0; j<m; ++j)
  {
    for(unsigned int i=0; i<n; ++i)
    {
      p[i] = path_action[i].eval(t[j]);
      dp[i] = path_action[i].deriv(t[j]);
    }
    lag_val[j] = lag(t[j], p, dp);
  }
  return lag_val;
}

template<typename T, typename Path, typename Lag> template<typename Gau_Kron>
void action<T, Path, Lag>::eval_helper(T left, T right, T D, T D_tol, T &integral, T &e)
{
  const Gau_Kron table;
  std::vector<T> fnodes(table.order, 0);
  T L1_lr=0; // L1 norm of local interval
  T D_lr;
  T scale_factor = (right - left)/2;
  T mid = (left+right)/2;
  T integral_l=0, integral_r=0; // left and right integral
  T e_l=0, e_r=0; // left and right error

  // calculate gauss and kronrod quadrature
  // calculate lagrangian at nodes
  fnodes[(table.order-1)/2] = mid;
  for(unsigned int i=0; i<(table.order-1)/2; ++i)
  {
    T tmp = (1.0+table.nodes[i])*scale_factor;
    fnodes[i] = right-tmp;
    fnodes[table.order-1-i] = tmp+left;
  }
  fnodes = eval_lagrangian(fnodes);

  // now integrate
  T int_kron=table.weight_kronrod[(table.order-1)/2]*fnodes[(table.order-1)/2];
  T int_gauss=0;
  // for gauss quadrature, zero only appears when order n = 4k+3
  if(table.order % 4 == 3)
  {
    int_gauss=table.weight_gauss[(table.order-3)/4]*fnodes[(table.order-1)/2];
    L1_lr = std::abs(int_gauss);
  }

  // gauss and kronrod quadrature formula
  for(unsigned int i=0; i<(table.order-1)/2; ++i)
  {
    int_kron += table.weight_kronrod[i]*(fnodes[i]+fnodes[table.order-1-i]);
    if(i %2 != 0)
    {
      int_gauss += table.weight_gauss[(i-1)/2]*\
      (fnodes[i]+fnodes[table.order-1-i]);
      L1_lr += table.weight_gauss[(i-1)/2]*\
      (std::abs(fnodes[i])+std::abs(fnodes[table.order-1-i]));
    }
  }

  // estimate error
  D_lr = int_kron - int_gauss;
  T length = std::abs(right-left);
  T tmp = std::abs(D_lr);
  T tmp_e = tmp*length;

  // integral is converged when
  // estimated relative error is lower than rtol
  if(tmp*std::sqrt(tmp*length) < D_tol*(eps2+L1_lr))
  {
    integral = scale_factor*int_kron;
    e = err_scale*tmp_e*std::sqrt(tmp_e); // update estimated error
    return;
  }

  // stop recurrsion due to tuncation
  // difference of Dlr and D is less than 
  // numerical epsilon 
  // (estimated by machine_eps*mean of abs value of function)
  if(std::abs(D_lr-D) < eps*(eps + L1_lr))
  {
    integral = scale_factor*int_kron;
    e = err_scale*tmp_e*std::sqrt(tmp_e);
    return;
  }

  // otherwise divide interval by half
  eval_helper<Gau_Kron>(left, mid, D_lr, D_tol, integral_l, e_l);
  eval_helper<Gau_Kron>(mid, right, D_lr, D_tol, integral_r, e_r);
  integral = integral_l + integral_r;
  e = e_l + e_r;
  return;
}

template<typename T, typename Path, typename Lag>
T action<T, Path, Lag>::eval_quadgk(T left, T right, unsigned int n, T &e)
{ 

  T D_tol = rtol/1000;
  T integral=0;

  // currently only supports n=15,21,31,41,51,61
  if(n==15)
  eval_helper<gau_kron_table<T, 15>>(left, right, 0, D_tol, integral, e);
  if(n==21)
  eval_helper<gau_kron_table<T, 21>>(left, right, 0, D_tol, integral, e);
  if(n==31)
  eval_helper<gau_kron_table<T, 31>>(left, right, 0, D_tol, integral, e);
  if(n==41)
  eval_helper<gau_kron_table<T, 41>>(left, right, 0, D_tol, integral, e);
  if(n==51)
  eval_helper<gau_kron_table<T, 51>>(left, right, 0, D_tol, integral, e);
  if(n==61)
  eval_helper<gau_kron_table<T, 61>>(left, right, 0, D_tol, integral, e);
  return integral;
}

template<typename T, typename Path, typename Lag>
T action<T, Path, Lag>::eval_qthsh(T left, T right, unsigned int max_order, T &e)
{
  T mid = (left+right)/2;
  T scale_coord = (right-left)/2;
  T scale_int = h_pi*scale_coord;
  T h = 2; // step_size
  T dt_pre = std::exp(T(1)); // previous step size
  T dt = dt_pre; // current step size
  T integral_pre = 0; // previous integration value
  T integral = eval_lagrangian(mid); // integration value
  T L1 = std::abs(integral); // L1 norm
  T tmp_e = 0; // error estimation

  for(unsigned int depth=0; depth<= max_order; ++depth)
  {
    T t = dt;
    T q = 0; 
    T f_l=0;
    T f_r=0;
    T corr = 0; // initialize correction term
    T L1_corr = 0; // correction term for L1
    h /= 2; // decrease step size by half
    do
    {
      T u = std::exp(h_pi*(1/t-t)); // exp(-pi*sinh(kh))
      T r = 2*u/(1+u); // (1-x_k)
      T dr = scale_coord*r;
      T w = (t+1/t)*r/(1+u);
      if(left+dr > left)
      f_l = eval_lagrangian(left+dr);
      if(right > right-dr)
      f_r = eval_lagrangian(right-dr);
      q = w*(f_l+f_r);
      corr += q;
      L1_corr += std::abs(q);
      t *= dt_pre;
    }
    // correction terms may not be changed
    // due to numerical tuncation
    while(std::abs(q) > eps*(eps+std::abs(corr)));
    integral_pre = integral; // previous term 
    integral += corr;
    tmp_e = std::abs(corr-integral_pre); //
    if(tmp_e < rtol*(eps2+L1))
    break; // integral meets tolerance
    if(std::abs(corr-integral_pre) < eps*(eps+std::abs(integral)))
    break; // numerical tuncation reaches
    dt_pre = dt;
    dt = std::sqrt(dt_pre);
    L1 += L1_corr; // update L1;
  }
  e = h*scale_int*tmp_e;
  return h*scale_int*integral;
}

template<typename T, typename Path, typename Lag>
T action<T, Path, Lag>::eval(T &e)
{
  if(! vaildity)
  return 0;
  T t0, t1;
  std::tie(t0, t1) = path_action[0].get_endtimes();
  e = 0;
  return eval_quadgk(t0, t1, 31, e);
}

template<typename T, typename Path, typename Lag>
T action<T, Path, Lag>::eval(int method, unsigned int n, T &e)
{
  if(! vaildity)
  return 0;
  T t0, t1;
  std::tie(t0, t1) = path_action[0].get_endtimes();
  e = 0;
  if(method == 0)
  return eval_quadgk(t0, t1, n, e);
  if(method == 1)
  return eval_qthsh(t0, t1, n, e);
  return 0; // if method is not equal to the one of 0 or 1.
}
}
