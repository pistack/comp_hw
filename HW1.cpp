/*
HW1
C++ code for homework1 in computer1 class in Yonsei University
Using explicit Euler Method to solve Kepler problem
*/

#include <cmath>
#include <tuple>
#include <vector>
#include <iostream>

using namespace std;

tuple<vector<double>, vector<double>> HW1(double t0, double t1,
					  int n, double y0, double y0p)
{ 

    // initialize the variable
    double spacing;
    double y_inv;
    vector<double> t(n+1, 0);
    vector<double> y(n+1, 0);

    spacing = (t1 - t0) / n;

    // use uniform n points bewteen t0 and t1
    // additional one point needed for end point
    for(int i=0; i< n+1; i++)
      {
        t[i] = i*spacing;
      }

    // initial condition
    y[0] = y0;
    y[1] = 1/n*(t1-t0)*y0p + y[0];

    // Solve 2nd order ODE using
    // Explict Euler Method

    for(int i = 2; i < n+1; i++)
      {
        y_inv = 1 / y[i-1];
        y[i] = pow(spacing, 2)*(pow(y_inv, 3) - pow(y_inv, 2)) + \
        2*y[i-1] - y[i-2];
      }

    return make_tuple(t, y);
}

int main(void) {
    double t0 = 0;
    double t1 = 10;
    int n = 100;
    double y0 = 0.9;
    double y0p = 0;
    vector<double> t, y;
    tie(t, y) = HW1(t0, t1, n, y0, y0p);

    for (int i = 0; i < n+1; i++)
      {
        cout << t[i] << '\t' << y[i] << endl;
      }

    return 0;
}
