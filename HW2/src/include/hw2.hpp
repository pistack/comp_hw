/*!
 * @file hw2.hpp
 * @brief headerfile for homework2 of Computer1 class in Yonsei University
 *        Use numerical integration to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 28.
 */

#ifndef HW2_H
#define HW2_H

#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>

 /// @brief HW2: Solve Kepler problem via numerical integration
 ///  from zeta_min to zeta_max
 /// If type T is not equal to one of float, double, long double
 /// then behavior of HW2 is undefined.
 /// @param zeta_min minimum value of zeta, 
 ///   for constraint motion 0.5 < zeta_min < 1
 /// @param t0 initial time
 /// @param n number of points to evaluate
 /// @return tuple of time and zeta
 /// @see \ref num_int
template <typename T>
std::tuple<std::vector<T>, std::vector<T>>
HW2(T zeta_min, T t0, int n);

#include "hw2.tpp"

#endif
