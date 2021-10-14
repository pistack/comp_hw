/*!
 * @file fourier_path.cpp
 * @brief evaluate path and derivative
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#include <cmath>
#include <algorithm>
#include <functional>
#include "fourier.hpp"

using namespace std;

fourier_path & fourier_path::operator=(const fourier_path & copy)
{
  // copy initial conditions
  p_t0 = copy.p_t0; p_tf = copy.p_tf; p_0 = copy.p_0; p_f = copy.p_f; 
  // copy fourier function, scaler and adder
	p_func = copy.p_func; p_vaild = copy.p_vaild; 
  scale = copy.scale; add = copy.add;
  return *this;
}

void fourier_path::init_helper()
{
  // init variable
  double tstp, tstf;
  vector<double> t_ends = {p_t0, p_tf};
  vector<double> f_ends(2, 0);

  f_ends = p_func.eval(t_ends);
  tstp = p_f - p_0;
  tstf = f_ends[1]-f_ends[0];
  if(abs(tstp)<1e-8 && abs(tstf)<1e-8)
  {
    p_vaild = true;
    scale = 1.0;
    add = p_0 - f_ends[0];
  }
  else if(abs(tstp)>1e-8 && abs(tstf)>1e-8)
  {
    p_vaild = true;
    scale = tstp/tstf;
    add = p_0 - scale*f_ends[0];
  }
  else
  {
    p_vaild = false;
    scale = 0;
    add = 0;
  }
}

void fourier_path::update(fourier fourier)
{
  p_func = fourier;
  init_helper();
}

bool fourier_path::is_vaild()
{
  return p_vaild;
}

double fourier_path::eval(double t)
{
  return scale*p_func.eval(t)+add;
}

vector<double> fourier_path::eval(vector<double> t)
{
  // init variable
  int n = t.size();
  vector<double> result(n, 0);
  
  result = p_func.eval(t);

  // scale first
  transform(result.begin(), result.end(), result.begin(),
	    bind2nd(multiplies<double>(), scale));

  // and then add
  transform(result.begin(), result.end(), result.begin(),
	    bind2nd(plus<double>(), add));

  return result;
}

double fourier_path::deriv(double t)
{
  return scale*p_func.deriv(t);
}

vector<double> fourier_path::deriv(vector<double> t)
{
  // init variable
  int n = t.size();
  vector<double> result(n, 0);
  
  result = p_func.deriv(t);

  transform(result.begin(), result.end(), result.begin(),
	    bind2nd(multiplies<double>(), scale));

  return result;
}

double fourier_path::nderiv(int n, double t)
{
  return scale*p_func.nderiv(n, t);
}

vector<double> fourier_path::nderiv(int n, vector<double> t)
{
  // init variable
  int n_t = t.size();
  vector<double> result(n_t, 0);
  
  result = p_func.nderiv(n, t);

  transform(result.begin(), result.end(), result.begin(),
	    bind2nd(multiplies<double>(), scale));

  return result;
}
  
