/*!
 * @file hw1.cpp
 * @brief code for homework1 of Computer1 class in Yonsei University
 * Use finite difference method to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 10.
 */

#include "hw1.hpp"

using namespace std;

tuple<vector<double>, vector<double>> HW1(double t0, double t1,
					  int n, double y0, double y0p)
{ 

    // initialize the variable
    double spacing;
    double y_inv;
    vector<double> t(n+1, 0);
    vector<double> y(n+1, 0);

    spacing = (t1 - t0) / double(n);

    // use uniform n points bewteen t0 and t1
    // additional one point needed for end point
    for(int i=0; i< n+1; i++)
      {
        t[i] = i*spacing;
      }

    // initial condition
    y[0] = y0;
    // estimate y1 using 2nd order talyor expension
    y[1] = spacing*(y0p+0.5*spacing*
    (pow(1.0/y0, 3.0) - pow(1.0/y0, 2.0))) + \
    y[0];

    // Solve 2nd order ODE using
    // Explict Euler Method

    for(int i = 2; i < n+1; i++)
      {
        y_inv = 1 / y[i-1];
        y[i] = pow(spacing, 2.0) * \
        (pow(y_inv, 3.0) - pow(y_inv, 2.0)) + \
        2*y[i-1] - y[i-2];
      }

    return make_tuple(t, y);
}
