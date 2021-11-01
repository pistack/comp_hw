/*!
 * @file action.tpp
 * @ingroup libfourier
 * @brief evaluates action
 * @author pistack (Junho Lee)
 * @date 2021. 11. 1.
 */

template<typename T, typename Lag>
action<T, Lag> & action<T, Lag>::operator=(const action<T, Lag> &copy)
{
  atol = copy.atol;
  path_action = copy.path_action;
  vaildity = copy.vaildity;
  return *this;
}

/// check vaildity 
template<typename T, typename Lag>
void action<T, Lag>::check_vaild()
{
  int n = path_action.size();
  vaildity = false;

  for(int i=0; i<n; i++)
  {
    if(! path_action[i].is_vaild())
    return;
  }
  vaildity = true;
  return;
}

/// update method 
template<typename T, typename Lag>
void action<T, Lag>::update(std::vector<fourier_path<T>> path)
{path_action = path; check_vaild();}

template<typename T, typename Lag>
void action<T, Lag>::update(T abs_tol)
{atol=abs_tol;}

template<typename T, typename Lag>
void action<T, Lag>::update(std::vector<fourier_path<T>> path, T abs_tol)
{atol=abs_tol; path_action = path; check_vaild();}

/// get vaildity 
template<typename T, typename Lag>
bool action<T, Lag>::is_vaild()
{return vaildity;}

template<typename T, typename Lag>
T action<T, Lag>::eval_lagrangian(T t)
{
  int n = path_action.size();
  std::vector<T> p(n, 0); // path
  std::vector<T> dp(n, 0); // derivative of path
  Lag lag; // lagrangian
  for(int i=0; i<n; i++)
  {
    p[i] = path_action[i].eval(t);
    dp[i] = path_action[i].deriv(t);
  }
  return lag(t, p, dp);
}

template<typename T, typename Lag>
std::vector<T> action<T, Lag>::eval_lagrangian(std::vector<T> t)
{
  int n = path_action.size();
  int m = t.size();
  std::vector<T> p(n, 0); // path
  std::vector<T> dp(n, 0); // derivative of path
  std::vector<T> lag_val(m, 0); // result
  Lag lag; // lagrangian
  for(int j=0; j<m; j++)
  {
    for(int i=0; i<n; i++)
    {
      p[i] = path_action[i].eval(t[j]);
      dp[i] = path_action[i].deriv(t[j]);
    }
    lag_val[j] = lag(t[j], p, dp);
  }
  return lag_val;
}

template<typename T, typename Lag> template<typename Gau_Kron>
T action<T, Lag>::eval_helper(T left, T right, T D, T D_tol)
{
  Gau_Kron table;
  std::vector<T> tnodes(table.order, 0);
  std::vector<T> fnodes(table.order, 0);
  T eps = 10.0*std::numeric_limits<T>::epsilon();
  T D_lr;
  T Delta;
  T scale_factor = 0.5*(right - left);
  T mid = 0.5*(left+right);
  tnodes[(table.order-1)/2] = mid;
  for(int i=0; i<(table.order-1)/2; i++)
  {
    T tmp = (1.0+table.nodes[i])*scale_factor;
    tnodes[i] = right-tmp;
    tnodes[table.order-1-i] = tmp+left;
  }
  fnodes = eval_lagrangian(tnodes);
  T int_kron=table.weight_kronrod[(table.order-1)/2]*fnodes[(table.order-1)/2];
  T int_gauss=table.weight_gauss[(table.order-3)/4]*fnodes[(table.order-1)/2];
  for(int i=0; i<(table.order-1)/2; i++)
  {
    int_kron += table.weight_kronrod[i]*(fnodes[i]+fnodes[table.order-1-i]);
    if(i %2 != 0)
    int_gauss += table.weight_gauss[(i-1)/2]*(fnodes[i]+fnodes[table.order-1-i]);
  }

  D_lr = int_kron - int_gauss;
  T tmp = std::abs(D_lr);
  if(tmp*std::sqrt(tmp*(right-left)) < D_tol)
  return scale_factor*int_kron;

  /// stop recurrsion due to tuncation
  /// 1. length of interval is less than
  /// numerical epsilon
  /// 2. difference of Dlr and D is less than 
  /// numerical epsilon

  if(right-left < eps*(eps+0.5*(std::abs(right)+std::abs(left))) ||
  std::abs(D_lr-D) <
  eps*(eps+std::abs(int_kron)))
  return scale_factor*int_kron;

  // otherwise divide interval by half
  return eval_helper<Gau_Kron>(left, mid, D_lr, D_tol) + \
  eval_helper<Gau_Kron>(mid, right, D_lr, D_tol);
}

template<typename T, typename Lag>
T action<T, Lag>::eval_quadgk(T left, T right, int n)
{ 

  T D_tol = atol/1000.0/(right-left);

  if(n==15)
  return eval_helper<gau_kron_table<T, 15>>(left, right, 0.0, D_tol);
  if(n==21)
  return eval_helper<gau_kron_table<T, 21>>(left, right, 0.0, D_tol);
  if(n==31)
  return eval_helper<gau_kron_table<T, 31>>(left, right, 0.0, D_tol);
  if(n==41)
  return eval_helper<gau_kron_table<T, 41>>(left, right, 0.0, D_tol);
  if(n==51)
  return eval_helper<gau_kron_table<T, 51>>(left, right, 0.0, D_tol);
  if(n==61)
  return eval_helper<gau_kron_table<T, 61>>(left, right, 0.0, D_tol);
  return 0; // un supported order
}

template<typename T, typename Lag>
T action<T, Lag>::eval()
{
  if(! vaildity)
  return 0;
  T t0, t1;
  std::tie(t0, t1) = path_action[0].get_endtimes();
  return eval_quadgk(t0, t1, 31);
}

