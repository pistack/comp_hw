/*!
 * @file math_const.hpp
 * @brief provides mathematical constent with different precision
 * @author pistack (Junho Lee)
 * @date 2021. 11. 15.
*/

#ifndef MATH_CONST_H
#define MATH_CONST_H

#include <limits>

namespace libpath {
/// @brief provides pi
/// @param T precision should be one of
/// float, double or long double 
template<typename T>
constexpr T PI();

/// @brief provides exponential constant
/// @param T precision should be one of
/// float, double or long double
template<typename T>
constexpr T EXP1();
}


#include "math_const/math_const.tpp"

#endif