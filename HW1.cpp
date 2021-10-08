/*
HW1
C++ code for homework1 in computer1 class in Yonsei University
Using explicit Euler Method to solve Kepler problem
*/

#include <tuple>
#include <vector>
#include <iostream>

using namespace std;

tuple<vector<double>, vector<double>> HW1(double t0, double t1,
int n, double y0, double y0p) { 

    // initialize the variable
    vector<double> t(n, 0);
    vector<double> y(n, 0);

    // use uniform n points bewteen t0 and t1
    // for time
    for(int i=0; i< n; i++) {
        t[i] = (t1-t0)/n*i;
    }

    // initial condition
    y[0] = y0;
    y[1] = 1/n*(t1-t0)*y0p + y[0];

    // Solve 2nd order ODE using
    // Explict Euler Method

    for(int i = 0; i < n-2; i++) {
        y[i+2] = ((t1-t0)/n)**2*(1/y[i]**3 - 1/y[i]**2) + \
        2*y[i+1] - y[i];
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

    for(int i = 0, i< n; i++) {
        cout << t[i] << ',' << y[i] << endl;
    }

    return 0;
}