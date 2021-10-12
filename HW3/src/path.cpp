/*!
 * @file path.cpp
 * @brief evaluate path and derivative
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#include "basic.hpp"

using namespace std;

void fourier_path::get_init(double t_init, double t_fin,
  double init, double fin, double period, vector<double> &c)
{
  // init variable
  vector<double> t_ends = {t_init, t_fin};
  vector<double> f_ends(2, 0);

  p_init = init; 
  p_final = fin;
  p_func.init(c.size()/2, period, c);
  f_ends = p_func.eval(t_ends);
  scale = (p_final - p_init)/(f_ends[1]-f_ends[0]);
  add = p_init - scale*f_ends[0];
}
