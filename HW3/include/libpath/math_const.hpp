/*!
 * @file math_const.hpp
 * @brief provides mathematical constent with different precision
 * @author pistack (Junho Lee)
 * @date 2021. 11. 15.
*/

#ifndef MATH_CONST_H
#define MATH_CONST_H

#include <limits>
#include <iostream>

namespace libpath {
/// @brief provides pi
/// @param T precision should be one of
/// float, double or long double 
template<typename T>
constexpr T PI();
}


#include "math_const/math_const.tpp"

#endif