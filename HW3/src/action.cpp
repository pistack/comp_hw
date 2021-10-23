/*!
 * @file action.cpp
 * @ingroup fourier
 * @brief evaluates action
 * @author pistack (Junho Lee)
 * @date 2021. 10. 23.
 */

#include <cmath>
#include <cerrno>
#include <vector>
#include "action.hpp"

using namespace std;

action & action::operator=(const action &copy)
{
  atol = copy.atol; rtol = copy.rtol;
  path_action = copy.path_action; lagranian = copy.lagranian;
  vaildity = copy.vaildity;
  return *this;
}

void action::check_vaild()
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

void action::update(vector<fourier_path> path)
{path_action = path; check_vaild();}

bool action::is_vaild()
{return vaildity;}

double action::eval_lagranian(double t)
{
  int n = path_action.size();
  vector<double> p(n, 0); // path
  vector<double> dp(n, 0); // derivative of path
  for(int i=0; i<n; i++)
  {
    p[i] = path_action[i].eval(t);
    dp[i] = path_action[i].deriv(t);
  }
  return lagranian(t, p, dp);
}

double action::eval_helper(double left, double mid, double right, 
double fleft, double fmid, double fright, 
double integral, double tol, int depth)
{
  double eps = 15*(tol + rtol*abs(integral));
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
  double lmid = 0.5*(left+mid);
  double rmid = 0.5*(right+mid);
  double flmid = eval_lagranian(lmid);
  double frmid = eval_lagranian(rmid);
  double stepsize_l = lmid - left;
  double stepsize_r = right - rmid;
  double integral_l = stepsize_l/3.0*(fleft+4.0*flmid+fmid);
  double integral_r = stepsize_r/3.0*(fright+4.0*frmid+fmid);
  double delta = integral_l + integral_r - integral;
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

double action::eval()
{
  if(! vaildity)
  {
    errno = EINVAL;
    return 0;
  }
  double t_0, t_1;
  tie(t_0, t_1) = path_action[0].get_endtimes();
  double tmid = 0.5*(t_0+t_1);
  double fleft = eval_lagranian(t_0);
  double fmid = eval_lagranian(tmid);
  double fright = eval_lagranian(t_1);
  double stepsize = tmid - t_0;
  double integral = stepsize/3.0*(fleft+4.0*fmid+fright);

  return eval_helper(t_0, tmid, t_1, fleft, fmid, fright,
  integral, atol, 0);
}

