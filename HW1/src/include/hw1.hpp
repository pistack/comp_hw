/*!
 * @file hw1.hpp
 * @brief Header file for homework1 of Computer1 class in Yonsei University
 * Use explicit Euler Method to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 10.
 */

#include <cmath>
#include <tuple>
#include <vector>


/*!
 * @brief HW1: Solve Kepler problem via explicit Euler Method
 *        with inital condition 
 *        - zeta(0) = z_0
 *        - zeta'(0) = z'_0
 *        see HW1.pdf for futher detail
 * @param t0 initial time
 * @param t1 final time
 * @param n number of gird points to evaluate
 * @param y0 initial condition for zeta(0)
 * @param y0p intial condition for zeta'(0)
 * @return tuple of time and zeta
 */
std::tuple<std::vector<double>, std::vector<double>>
HW1(double t0, double t1, int n, double y0, double y0p);

