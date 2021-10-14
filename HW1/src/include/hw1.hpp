/*!
 * @file hw1.hpp
 * @brief Header file for homework1 of Computer1 class in Yonsei University
 * Use finite difference method to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 10.
 */

#ifndef HW1_H
#define HW1_H
#endif

#include <cmath>
#include <tuple>
#include <vector>

 /// @brief HW1: Solve Kepler problem via finite difference Method
 /// @param t0 initial time
 /// @param t1 final time
 /// @param n number of gird points to evaluate
 /// @param y0 initial condition for zeta
 /// @param y0p intial condition for derivative of zeta
 /// @param theta0 initial condition for theta
 /// @return tuple of time, zeta and theta
 /// @see \ref fdm
 /// @see \ref theta
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
HW1(double t0, double t1, int n, double y0, double y0p, double theta0);

