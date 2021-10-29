/*!
 * @file action.cpp
 * @ingroup libfourier
 * @brief evaluates action<T>
 * @author pistack (Junho Lee)
 * @date 2021. 10. 29.
 */

template<typename T, typename Lag>
action<T, Lag> & action<T, Lag>::operator=(const action<T, Lag> &copy)
{
  atol = copy.atol; rtol = copy.rtol;
  path_action = copy.path_action;
  vaildity = copy.vaildity;
  return *this;
}

template<typename T, typename Lag>
void action<T, Lag>::check_vaild()
{
  if(!path_action[0].is_vaild())
  return;

  int n = path_action.size();

  for(int i=1; i<n; i++)
  {
    if(! path_action[i].is_vaild() ||
    (path_action[i-1].get_endtimes() !=
    path_action[i].get_endtimes()))
    return;
  }
  vaildity = true;
  return;
}

template<typename T, typename Lag>
void action<T, Lag>::update(std::vector<fourier_path<T>> path)
{path_action = path; check_vaild();}

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

template<typename T, typename Lag>
T action<T, Lag>::eval_helper(T left, T mid, T right, 
T fleft, T fmid, T fright, 
T integral, T tol, int depth)
{
  T eps = 15*(tol + rtol*abs(integral));
  if(depth > MAXDEPTH)
  {
    errno = ERANGE;
    return integral;
  }
  if(eps == eps/2 || left == mid)
  {
    errno = EDOM;
    return integral;
  }
  T lmid = 0.5*(left+mid);
  T rmid = 0.5*(right+mid);
  T flmid = eval_lagranian(lmid);
  T frmid = eval_lagranian(rmid);
  T stepsize_l = lmid - left;
  T stepsize_r = right - rmid;
  T integral_l = stepsize_l/3.0*(fleft+4.0*flmid+fmid);
  T integral_r = stepsize_r/3.0*(fright+4.0*frmid+fmid);
  T delta = integral_l + integral_r - integral;
  if(abs(delta) > eps || depth==0)
  {
    integral_l = eval_helper(left, lmid, mid,
    fleft, flmid, fmid, integral_l, tol/2, depth+1);
    integral_r = eval_helper(mid, rmid, right,
    fmid, frmid, fright, integral_r, tol/2, depth+1);
    return integral_l + integral_r;
  }
  else
  return integral_l + integral_r + delta/15;
}

template<typename T, typename Lag>
T action<T, Lag>::eval()
{
  if(! vaildity)
  {
    errno = EINVAL;
    return 0;
  }
  T t_0, t_1;
  std::tie(t_0, t_1) = path_action[0].get_endtimes();
  T tmid = 0.5*(t_0+t_1);
  T fleft = eval_lagranian(t_0);
  T fmid = eval_lagranian(tmid);
  T fright = eval_lagranian(t_1);
  T stepsize = tmid - t_0;
  T integral = stepsize/3.0*(fleft+4.0*fmid+fright);

  return eval_helper(t_0, tmid, t_1, fleft, fmid, fright,
  integral, atol, 0);
}

