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
  for(int i=0; i<n; i++)
  {
    p[i] = path_action[i].eval(t);
    dp[i] = path_action[i].deriv(t);
  }
  Lag lag; // lagrangian
  return lag(t, p, dp);
}

/// evaluation helper
template<typename T, typename Lag>
T action<T, Lag>::eval_helper(T left, T mid, T right, 
T fleft, T fmid, T fright, 
T integral, T D)
{
  T lmid = 0.5*(left+mid);
  T rmid = 0.5*(right+mid);
  T flmid = eval_lagranian(lmid);
  T frmid = eval_lagranian(rmid);
  T stepsize_l = lmid - left;
  T stepsize_r = right - rmid;
  T integral_l = stepsize_l/3.0*(fleft+4.0*flmid+fmid);
  T integral_r = stepsize_r/3.0*(fright+4.0*frmid+fmid);
  T delta = integral_l + integral_r - integral;
  T Dlr = 0.25*(fleft - 4.0*flmid + 6.0*fmid - 4.0*frmid + fright);
  T absD = std::abs(D);
  T absDlr = std::abs(Dlr);

  /// if interval is converged than
  /// 15 epsilon acceleration is used 
  if(absDlr < D_tol && absDlr < absD)
  return integral_l + integral_r + delta/15;

  /// turn off 15 epsilon acceleration
  /// when difference of Dlr and D is less than 
  /// numerical tuncation
  if(std::abs(Dlr-D) <
  50.0*std::numeric_limits<T>::epsilon()*(absD+absDlr))
  return integral_l + integral_r;

  // otherwise, divide interval by two subintervals.
  integral_l = eval_helper(left, lmid, mid, fleft, flmid, fmid,
  integral_l, Dlr);
  integral_r = eval_helper(mid, rmid, right, fmid, frmid, fright,
  integral_r, Dlr);

  return integral_l + integral_r;
}

template<typename T, typename Lag>
T action<T, Lag>::eval()
{
  if(! vaildity)
  {
    errno = EINVAL;
    return 0;
  }
  T t0, t1;
  std::tie(t0, t1) = path_action[0].get_endtimes();
  D_tol = 180.0*atol/(t1-t0);
  T tmid = 0.5*(t0+t1);
  T fleft = eval_lagranian(t0);
  T fmid = eval_lagranian(tmid);
  T fright = eval_lagranian(t1);
  T stepsize = tmid - t0;
  T integral = stepsize/3.0*(fleft+4.0*fmid+fright);

  return eval_helper(t0, tmid, t1, fleft, fmid, fright,
  integral, 0.0);
}

