/*!
 * @file bezier_path.hpp
 * @ingroup libpath
 * @brief headerfile for path approximated by bezier curve
 * @author pistack (Junho Lee)
 * @date 2021. 11. 10.
 */

#ifndef BEZIER_PATH_H
#define BEZIER_PATH_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <vector>
#include "bezier.hpp"

namespace libpath {

/// @ingroup libpath
/// @brief Class which approximate path by \f$ n \f$ th order Bezier curve.
/// \f{equation}{\label{eq:bezier_path}
/// \Phi(t) = \mathrm{scale}_2 \cdot 
/// B\left(\frac{t-t_i}{t_f-t_i}, n, 
/// \{c_0, c_1, \dotsc, \mathrm{scale}_1 \cdot c_n\}\right)
/// \f}
/// ,where \f$ c_i \f$ is control points of bezier curve,
/// \f$ \mathrm{scale}_1 \f$ is the first scaling parameter to match ratio of
/// initial and final value of path and \f$ \mathrm{scale}_2 \f$ is
/// the second scaling parameter to match boundary condition.
/// @param T precision should be one of float, double or long double.
/// @note If initial value of path, given by boundary condition, is zero
/// then it sets \f$ c_0 = 0 \f$.
/// @warning If bezier curve is not appropriate to approximate path then
/// eval and deriv method returns aways zero.
template<typename T>
class bezier_path
{
    private:
    // initial condition
    T t_init, t_final, p_init, p_final;
    bezier<T> B; // bezier curve
    bool vaildity; // vaildity of path
    T scale2; // second scaling parameter

    /// @brief helper function for 
    /// initialization of bezier_path class
    void init_helper();

    public:
    /// @brief initialize bezier_path class
    bezier_path() {}

    /// @brief initialize bezier_path class
    /// @param t_init_ initial time
    /// @param t_final_ final time
    /// @param p_init_ initial value of path
    /// @param p_final_ final value of path
    bezier_path(T t_init_, T t_final_, T p_init_, T p_final_)
    : t_init(t_init_), t_final(t_final_), p_init(p_init_), p_final(p_final_)
    {}

    /// @brief initialize bezier_path class
    /// @param t_init_ initial time
    /// @param t_final_ final time
    /// @param p_init_ initial value of path
    /// @param p_final_ final value of path
    /// @param B_ Bezier curve used to approximate path
    bezier_path(T t_init_, T t_final_, T p_init_, T p_final_, bezier<T> B_)
    : t_init(t_init_), t_final(t_final_), p_init(p_init_), p_final(p_final_), \
    B(B_)
    {init_helper();}

    /// @brief copy constructor of bezier_path class
    bezier_path(const bezier_path<T> & copy)
    : t_init(copy.t_init), t_final(copy.t_final), \
    p_init(copy.p_init), p_final(copy.p_final), B(copy.B), \
    vaildity(copy.vaildity), scale2(copy.scale2)
    {}

    /// @brief overloading of assignment operator for 
	/// bezier_path class
	bezier_path<T> & operator=(const bezier_path<T> &copy)
    {
        t_init = copy.t_init; t_final = copy.t_final;
        p_init = copy.p_init; p_final = copy.p_final;
        B = copy.B; vaildity = copy.vaildity;
        scale2 = copy.scale2;
        return *this;
    }

    /// @brief update bezier curve
	/// @param B_ bezier curve to update
	void update(bezier<T> B_)
    {B = B_; init_helper();}

    /// @brief check whether or not the bezier curve is valid to
	/// approximate path
	/// @return vaildity of bezier curve
	bool is_vaild() const
    {return vaildity;}

	/// @brief get initial and final time of path
	/// @return tuple of initial and final time of path
	std::tuple<T, T> get_endtimes() const
    {return std::make_tuple(t_init, t_final);}

	/// @brief get modified control points
	/// @return modified control points
	std::vector<T> get_ctrl_pts() const
    {return B.get_crtl_pts();}

	/// @brief get the second scaler parameter
	/// @return second scaler parameter
	T get_scaler() const
    {return scale2;}

    /// @brief evaluate path approximated by bezier curve at given point
    /// @param t point at which bezier_path is evaluated
    T eval(T t)
    {return scale2*B.eval((t-t_init)/(t_final-t_init));}

    /// @brief evaluate path approximated by bezier curve at given points
    /// @param t points at which bezier_path is evaluated
    std::vector<T> eval(std::vector<T> t);

    /// @brief evaluate derivative of bezier_path at given point
    /// @param t point at which bezier_path is evaluated
    T deriv(T t)
    {return scale2*B.deriv((t-t_init)/(t_final-t_init))/(t_final-t_init);}

    /// @brief evaluate derivative of bezier_path at given points
    /// @param t points at which bezier_path is evaluated
    std::vector<T> deriv(std::vector<T> t);
};

}

#include "bezier_path/bezier_path.tpp"

#endif