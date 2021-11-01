/*!
 * @file action.tpp
 * @ingroup libfourier
 * @brief evaluates action
 * @author pistack (Junho Lee)
 * @date 2021. 10. 31.
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

template<typename T, typename Lag>
T action<T, Lag>::eval_helper(T left, T right, T D)
{
  std::vector<T> tnodes(31, 0);
  std::vector<T> fnodes(31, 0);
  T D_lr;
  T Delta;
  T scale_factor = 0.5*(right - left);
  T mid = 0.5*(left+right);
  tnodes[15] = mid;
  for(int i=0; i<15; i++)
  {
    T tmp = (1.0+nodes[i])*scale_factor;
    tnodes[i] = right-tmp;
    tnodes[30-i] = tmp+left;
  }
  fnodes = eval_lagrangian(tnodes);
  T int_kron=weight_kronrod[15]*fnodes[15];
  T int_gauss=weight_gauss[7]*fnodes[15];
  for(int i=0; i<15; i++)
  {
    int_kron += weight_kronrod[i]*(fnodes[i]+fnodes[30-i]);
    if(i %2 != 0)
    int_gauss += weight_gauss[(i-1)/2]*(fnodes[i]+fnodes[30-i]);
  }

  D_lr = int_kron - int_gauss;
  T tmp = std::abs(D_lr);
  if(tmp*std::sqrt(tmp*(right-left)) < D_tol)
  return scale_factor*int_kron;

  /// when difference of Dlr and D is less than 
  /// numerical epsilon stop recurrsion
  if(std::abs(D_lr-D) <
  1.0*std::numeric_limits<T>::epsilon()*
  (std::numeric_limits<T>::min()+std::abs(int_kron)))
  return scale_factor*int_kron;

  // otherwise divide interval by half
  return eval_helper(left, mid, D_lr) + \
  eval_helper(mid, right, D_lr);
}

template<typename T, typename Lag>
T action<T, Lag>::eval_quadgk(T left, T right)
{ 

  D_tol = atol/1000.0/(right-left);

  return eval_helper(left, right, 0.0);
}

template<typename T, typename Lag>
T action<T, Lag>::eval()
{
  if(! vaildity)
  return 0;
  T t0, t1;
  std::tie(t0, t1) = path_action[0].get_endtimes();
  return eval_quadgk(t0, t1);
}

