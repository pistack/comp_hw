/*!
 * @file path.cpp
 * @brief evaluate path and derivative
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#include "basic.hpp"

using namespace std;

tuple<vector<double>, vector<double>>
eval_path(vector<double> &t, double init, double fin,
	  vector<double> &c)
{
  // initialize variable
  int n = t.size();
  int num_fourier = c.size()/2; // one term = one sine + one cosine
  double scale;
  double add;
  vector<double> path(n, 0);
  vector<double> deriv_path(n, 0);

  tie(path, deriv_path) = sum_of_fourier(t, c, num_fourier);
  scale = (fin-init)/(path[n-1]-path[0]);
  add = init - scale*path[0];
  path = scale_and_add_vector(path, scale, add);
  deriv_path = scale_and_add_vector(deriv_path, scale, 0.0);

  return make_tuple(path, deriv_path);
}
