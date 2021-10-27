/*!
 * @file hw1.hpp
 * @brief Header file for homework1 of Computer1 class in Yonsei University
 * Use finite difference method to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 28.
 */

#ifndef HW1_H
#define HW1_H

#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>

 /// @brief HW1: Solve Kepler problem via finite difference Method
 /// Behavior of HW1 is undefined when type T is not equal to one of
 /// float, double, long double.
 /// @param t0 initial time
 /// @param t1 final time
 /// @param n number of gird points to evaluate
 /// @param y0 initial condition for zeta
 /// @param y0p intial condition for derivative of zeta
 /// @param theta0 initial condition for theta
 /// @return tuple of time, zeta and theta
 /// @see \ref fdm
 /// @see \ref theta
template <typename T>
std::tuple<std::vector<T>, std::vector<T>, std::vector<T>>
HW1(T t0, T t1, int n, T y0, T y0p, T theta0);

#include "hw1.tpp"

#endif

