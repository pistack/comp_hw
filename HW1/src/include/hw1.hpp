/*!
 * @file hw1.hpp
 * @brief Header file for homework1 of Computer1 class in Yonsei University
 * Use explicit Euler Method to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 10.
 */

#ifndef HW1_H
#define HW1_H
#endif

#include <cmath>
#include <tuple>
#include <vector>

////////////////////////////////////////////////////////////////////

/*!
 * @mainpage Computer Homework 1
 * Solve Kepler problem via explicti Euler Method
 * @section req Requirements
 * To install this program, you should have
 * - C++ compiler like g++
 * - gnu make
 * @section install Installation
 * Type make, then we can see hw1 executable file in bin directory
 * @section hwoto How To Use
 * Execute hw1 then, it will interactively read
 * - inital condition
 * - number of gird points to evaluate
 * - output file name
 *
 * Then it computes and saves solution to file.
 * You can plot the result using usual plotting software like gnuplot
 */

///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////

///@page page1 Mathematics
///@tableofcontents
///@section math Mathematics
/// To solve Kepler problem, we need to solve
/// \f{equation}{\label{eq:ode_2nd}
/// \frac{\mathrm{d}^2 \zeta}{\mathrm{d}t^2} =
/// \frac{1}{\zeta^3} - \frac{1}{\zeta^2}
/// \f}
/// with initial condition
/// \f{align*}{
/// \zeta(t_0) &= \zeta_0 \\
/// \zeta'(t_0) &= \zeta'_0
/// \f}
/// To solve above 2nd order ordinary differential equation
/// \latexonly \eqref{eq:ode_2nd} \endlatexonly
/// we need to aproximate 

///////////////////////////////////////////////////////////////////////

/*!
 * @brief HW1: Solve Kepler problem via explicit Euler Method
 * @param t0 initial time
 * @param t1 final time
 * @param n number of gird points to evaluate
 * @param y0 initial condition for zeta
 * @param y0p intial condition for derivative of zeta
 * @return tuple of time and zeta
 */

std::tuple<std::vector<double>, std::vector<double>>
HW1(double t0, double t1, int n, double y0, double y0p);

