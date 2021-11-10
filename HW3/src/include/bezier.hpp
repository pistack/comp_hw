/*!
 * @file bezier.hpp
 * @ingroup libpath
 * @brief headerfile for bezier curve and path approximated by bezier curve
 * @author pistack (Junho Lee)
 * @date 2021. 11. 10.
 */

#ifndef BEZIER_H
#define BEZIER_H

#include <cstddef>
#include <vector>

namespace libpath {
/// @ingroup libpath
/// @brief Class which defines \f$ n \f$ th order Bezier curve.
/// \f{equation}{\label{eq:bezier}
/// B(t) = \sum_{i=0}^n c_i \binom{n}{i} (1-t)^{n-i}t^i
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
    bezier() {}

    /// @brief initialize bezier class
    /// @param n_ order of Bezier curve
    /// @param c_ control points of Bezier curve
    bezier(int n_, std::vector<T> c_)
    : n(n_), c(c_)
    {}

    /// @brief copy constructor of bezier class
    bezier(const bezier<T> & copy)
    : n(copy.n), c(copy.c)
    {}

    /// @brief overloading of assignment operator for 
	/// bezier class
	bezier<T> & operator=(const bezier<T> &copy)
    {
    n = copy.n; c = copy.c;
    return *this;
    }

    /// @brief get first control point
	/// @return first control point
	T get_first() const
    {return c[0];}

    /// @brief get last control point
	/// @return last control point
	T get_last() const
    {return c[n];}

    /// @brief get control points of bezier curve
    /// @return control points of bezier curve
    std::vector<T> get_crtl_pts() const
    {return c;}

    /// @brief update control points
	/// @param c_ control points to update
	void update(std::vector<T> c_)
    {c = c_;}

    /// @brief update first control point
	/// @param c_0 first control point to update
	void update_first(T c_0)
    {c[0] = c_0;}

    /// @brief update last control point
	/// @param c_n last control point to update
	void update_last(T c_n)
    {c[n] = c_n;}

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
}

#include "bezier/bezier.tpp"

#endif