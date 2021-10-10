/*!
 * @file hw2.hpp
 * @brief headerfile for homework2 of Computer1 class in Yonsei University
 *        Use numerical integration to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 10.
 */

#ifndef HW2_H
#define HW2_H
#endif

#include <cmath>
#include <tuple>
#include <vector>

/*!
 * @brief HW2: Solve Kepler problem via numerical integration
 *  from zeta_min to zeta_max
 *  see HW2.pdf for further detail.
 * @param zeta_min minimum value of zeta, 
    for constraint motion 0.5 < zeta_min < 1
 * @param t0 initial time
 * @param n number of points to evaluate
 * @return tuple of time and zeta
 */

std::tuple<std::vector<double>, std::vector<double>>
HW2(double zeta_min, double t0, int n);
