/*!
 * @file mcm_bezier.hpp
 * @ingroup libmcm
 * @brief headerfile for mcm_bezier class
 * which is derivated by libmcm::mcm class
 * @author pistack (Junho Lee)
 * @date 2021. 11. 12.
 */

#include "libpath/bezier_path.hpp"
#include "mcm.hpp"

namespace libmcm {
/// @brief Derivated class of libmcm::mcm class for
/// path approximated by bezier curve
/// @param T precision should be one of float, double, or
/// long double
/// @param Lag lagrangian 
/// functor class which has
/// time, path and derivative of path
/// as variable and it returns
/// value of lagranian at given time
/// @see libpath::bezier_path class
/// @see libmcm::mcm class
template<typename T, typename Lag>
class mcm_bezier :public mcm<T, libpath::bezier<T>, 
libpath::bezier_path<T>, Lag>
{
    public:
    /// @brief initalize mcm_bezier class
    mcm_bezier() {}

    /// @brief initalize mcm_bezier class
    /// @param t_0 initial time
    /// @param t_1 final time 
    /// @param p_0 value of path at initial time
    /// @param p_1 vale of path at final time
    /// @param abs_tol absolute tolerance for action integral
    /// @param n_ order of bezier curve
    mcm_bezier(T t_0, T t_1, std::vector<T> p_0, std::vector<T> p_1, 
    T abs_tol, unsigned int n_)
    : mcm<T, libpath::bezier<T>, libpath::bezier_path<T>,
    Lag>(t_0, t_1, p_0, p_1, abs_tol, n_)
    {}

    /// @brief copy constructor of mcm_bezier class
    mcm_bezier(const mcm_bezier<T, Lag> &copy)
    : mcm<T, libpath::bezier<T>, libpath::bezier_path<T>,
    Lag>(copy)
    {}

    /// @brief overloading of assignment operator for 
    /// mcm_bezier class
    mcm_bezier<T, Lag> & operator=(const mcm_bezier<T, Lag> &copy)
    {
        this -> mcm<T, libpath::bezier<T>, libpath::bezier_path<T>, 
        Lag>::operator=(copy);
        return *this;
    }

    /// @brief get scaler and coefficients of
    /// initial guess
    /// @return tuple of scaler and coefficients
    /// of initial guess
    std::tuple<std::vector<T>, std::vector<std::vector<T>>>
    get_init_coeff();

    /// @brief get scaler and coefficients of
    /// minimal guess
    /// @return tuple of scaler and coefficients
    /// of minimal guess
    std::tuple<std::vector<T>, std::vector<std::vector<T>>>
    get_min_coeff();

};
}
#include "mcm_bezier/mcm_bezier.tpp"