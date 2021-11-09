/*!
 * @file bezier.hpp
 * @ingroup libfourier
 * @brief headerfile for bezier curve
 * @author pistack (Junho Lee)
 * @date 2021. 11. 9.
 */

#ifndef BEZIER_H
#define BEZIER_H

#include <algorithm>
#include <vector>

/// @brief Class which defines \f$ n \f$ th order Bezier curve.
/// \f{equation}{\label{eq:bezier}
/// P(t) = \sum_{i=0}^n c_i \binom{n}{i} (1-t)^{n-i}t^i
/// \f}
/// ,where \f$ c_i \f$ is control points.
/// @param T precision should be one of float, double or long double.
/// @note For \f$ n \f$ th order Bezier curve, you should give
/// \f$ n+1 \f$ control points.
template<typename T>
class bezier
{
    private:
    int n; // order of Bezier curve
    std::vector<T> c; // control points of Bezier curve

    public:
    /// @brief initialize bezier class
    bezier() {};

    /// @brief initialize bezier class
    /// @param n_ order of Bezier curve
    /// @param c_ control points of Bezier curve
    bezier(int n_, std::vector<T> c_)
    : n(n_), c(c_)
    {};

    /// @brief copy constructor of bezier class
    bezier(const bezier<T> & copy)
    : n(copy.n), c(copy.c)
    {};

    /// @brief overloading of assignment operator for 
	/// bezier class
	bezier<T> & operator=(const bezier<T> &copy);

    /// @brief evaluate bezier curve at given point
    /// \f$ 0 \leq t \leq 1 \f$.
    /// @param t point at which bezier curve is evaluated
    T eval(T t);

    /// @brief evaluate bezier curve at given points
    /// \f$ 0 \leq t \leq 1 \f$.
    /// @param t points at which bezier curve is evaluated
    std::vector<T> eval(std::vector<T> t);

    /// @brief evaluate derivative of bezier curve at given point
    /// \f$ 0 \leq t \leq 1 \f$.
    /// @param t point at which bezier curve is evaluated
    T deriv(T t);

    /// @brief evaluate bezier curve at given points
    /// \f$ 0 \leq t \leq 1 \f$.
    /// @param t points at which bezier curve is evaluated
    std::vector<T> deriv(std::vector<T> t);
};

#include "bezier/bezier.tpp"

#endif