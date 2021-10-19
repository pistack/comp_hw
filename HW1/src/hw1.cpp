/*!
 * @file hw1.cpp
 * @brief code for homework1 of Computer1 class in Yonsei University
 * Use finite difference method to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 16.
 */

#include <algorithm>
#include "hw1.hpp"

using namespace std;

tuple<vector<double>, vector<double>, vector<double>> 
HW1(double t0, double t1,
int n, double y0, double y0p, double theta0)
{ 

    // initialize the variable
    double spacing;
    double y_inv;
    vector<double> t(n+1, 0);
    vector<double> y(n+1, 0);
    vector<double> theta(n+1, 0);

    spacing = (t1 - t0) / double(n);

    // use uniform n points bewteen t0 and t1
    // additional one point needed for end point
    // to avoid numerical problem add t0 later
    for(int i=1; i< n+1; i++)
      t[i] = t[i-1] + spacing;

    // initial condition
    y[0] = y0;
    theta[0] = theta0;
    // estimate y1 using 2nd order talyor expension
    y_inv = 1/y[0];
    spacing = t[1] - t[0];
    y[1] = spacing/y[0]*(y0p+0.5*spacing/y[0]*
    (y_inv - 1.0)) + y[0];
    // estimate theta1 using trapezoid rule
    theta[1] = theta[0] + \
    0.5*spacing*pow(y_inv, 2.0) * \
    (1.0+pow(y[0]/y[1], 2.0));

    // Solve 2nd order ODE using
    // Explict Euler Method

    for(int i = 2; i < n+1; i++)
    {
      spacing = (t[i] - t[i-2])/2.0; /// spacing for zeta
      y_inv = 1 / y[i-1];
      y[i] = pow(spacing/y[i-1], 2.0) * \
      (y_inv - 1.0) + \
      (2*y[i-1] - y[i-2]);
      spacing = t[i] - t[i-1]; // spacing for theta
      theta[i] = theta[i-1] + \
      0.5*spacing*pow(y_inv, 2.0) * \
      (1.0+pow(y[i-1]/y[i], 2.0));
    }

    transform(t.begin(), t.end(), t.begin(),
    [t0](double &x){return x += t0;});
    t[n] = t1;
    
    return make_tuple(t, y, theta);
}
