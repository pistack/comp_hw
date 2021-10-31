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
T action<T, Lag>::eval_lagranian(T t)
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
std::vector<T> action<T, Lag>::eval_lagranian(std::vector<T> t)
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
  std::vector<T> nodes{
    0.991455371120813, 0.949107912342759, 0.864864423359769,
    0.741531185599394, 0.586087235467691, 0.405845151377397, 
    0.207784955007898
  };
  std::vector<T> weight_gauss{
    0.129484966168870, 0.279705391489277, 0.381830050505119,
    0.417959183673469
  };
  std::vector<T> weight_kronrod{
    0.022935322010529, 0.063092092629979, 0.104790010322250,
    0.140653259715525, 0.169004726639267, 0.190350578064785,
    0.204432940075298, 0.209482141084728 
  };
  std::vector<T> tnodes(15, 0);
  std::vector<T> fnodes(15, 0);
  T D_lr;
  T Delta;
  T scale_factor = 0.5*(right - left);
  tnodes[7] = 0.5*(left+right);
  for(int i=0; i<7; i++)
  {
    T tmp = (1.0+nodes[i])*scale_factor;
    tnodes[i] = right-tmp;
    tnodes[14-i] = tmp+left;
  }
  fnodes = eval_lagranian(tnodes);
  T int_kron=weight_kronrod[7]*fnodes[7];
  T int_gauss=weight_gauss[3]*fnodes[7];
  for(int i=0; i<7; i++)
  {
    int_kron += weight_kronrod[i]*(fnodes[i]+fnodes[14-i]);
    if(i %2 != 0)
    int_gauss += weight_gauss[(i-1)/2]*(fnodes[i]+fnodes[14-i]);
  }

  D_lr = int_kron - int_gauss;
  if(std::abs(D_lr) < D_tol)
  return scale_factor*int_kron;

  /// when difference of Dlr and D is less than 
  /// numerical epsilon stop interation
  if(std::abs(D_lr-D) <
  5.0*std::numeric_limits<T>::epsilon()*(2.0+std::abs(D)+std::abs(D_lr)))
  return scale_factor*int_kron;

  // otherwise divide interval by half
  return eval_helper(left, tnodes[7], D_lr) + \
  eval_helper(tnodes[7], right, D_lr);
}

template<typename T, typename Lag>
T action<T, Lag>::eval()
{
  if(! vaildity)
  return 0;
  
  T t0, t1;
  std::tie(t0, t1) = path_action[0].get_endtimes();
  D_tol = std::pow(atol, 1.0/1.5)/100.0/(t1-t0);

  return eval_helper(t0, t1, 0.0);
}

