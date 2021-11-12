/*!
 * @file mcm_fourier.hpp
 * @ingroup libmcm
 * @brief headerfile for mcm_fourier class
 * which is derivated by libmcm::mcm class
 * @author pistack (Junho Lee)
 * @date 2021. 11. 12.
 */

#include "libpath/fourier_path.hpp"
#include "mcm.hpp"

namespace libmcm {
/// @brief Derivated class of libmcm::mcm class for
/// path approximated by fourier function
/// @param T precision should be one of float, double, or
/// long double
/// @param Lag lagrangian 
/// functor class which has
/// time, path and derivative of path
/// as variable and it returns
/// value of lagranian at given time
/// @see libpath::fourier_path class
/// @see libmcm::mcm class
template<typename T, typename Lag>
class mcm_fourier :public mcm<T, libpath::fourier<T>, 
libpath::fourier_path<T>, Lag>
{
    public:
    /// @brief initalize mcm_fourier class
    mcm_fourier() {}

    /// @brief initalize mcm_fourier class
    /// @param t_0 initial time
    /// @param t_1 final time 
    /// @param p_0 value of path at initial time
    /// @param p_1 vale of path at final time
    /// @param abs_tol absolute tolerance for action integral
    /// @param n_f number of fourier function to use
    /// @param period period of fourier function
    mcm_fourier(T t_0, T t_1, std::vector<T> p_0, std::vector<T> p_1, 
    T abs_tol, unsigned int n_f, T period)
    : mcm<T, libpath::fourier<T>, libpath::fourier_path<T>,
    Lag>(t_0, t_1, p_0, p_1, abs_tol, n_f, period)
    {}

    /// @brief copy constructor of mcm_fourier class
    mcm_fourier(const mcm_fourier<T, Lag> &copy)
    : mcm<T, libpath::fourier<T>, libpath::fourier_path<T>,
    Lag>(copy)
    {}

    /// @brief overloading of assignment operator for 
    /// mcm_fourier class
    mcm_fourier<T, Lag> & operator=(const mcm_fourier<T, Lag> &copy)
    {
        this -> mcm<T, libpath::fourier<T>, libpath::fourier_path<T>, 
        Lag>::operator=(copy);
        return *this;
    }

    /// @brief get adder, scaler and coefficients of
    /// initial guess
    /// @return tuple of adder, scaler and coefficients
    /// of initial guess
    std::tuple<std::vector<T>, std::vector<T>,
    std::vector<std::vector<T>>>
    get_init_coeff();

    /// @brief get adder, scaler and coefficients of
    /// minimal guess
    /// @return tuple of adder, scaler and coefficients
    /// of minimal guess
    std::tuple<std::vector<T>, std::vector<T>,
    std::vector<std::vector<T>>>
    get_min_coeff();

};
}