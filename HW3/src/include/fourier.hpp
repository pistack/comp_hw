/*!
 * @file fourier.hpp
 * @ingroup libpath
 * @brief headerfile for fourier function
 * @author pistack (Junho Lee)
 * @date 2021. 11. 10.
 */

#ifndef FOURIER_H
#define FOURIER_H

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <vector>

namespace libpath {

template<typename T>
constexpr T pi = std::acos(T(-1)); ///< pi

/// @brief Class which defines fourier function
/// \f{equation}{\label{eq:fourier_func} 
/// \phi\{c, T\}(t) = \sum_{i=0}^{n_f-1} c_{2i} \sin\left(\frac{2i\pi}{T} t \right) +
/// c_{2i+1} \cos\left(\frac{2i\pi}{T} t \right)
/// \f}
/// , where \f$ n_f \f$ is the number of sine and cosine function used in
/// fourier function, \f$ c \f$ is the weight and \f$ T \f$ is the period of
/// fourier function.
/// @param T precision should be one of float, double,
/// and long double
/// @ingroup libpath

template<typename T>
class fourier
{
	private:

	int f_num_fourier; // number of sine and cosine function to add
	T f_period; // the period of fourier function 
	std::vector<T> f_c; // coefficients

	public:

	/// @brief initialize fourier class
	fourier() {}

	/// @brief initialize fourier class
	/// @param num_fourier number of sine and consine function to add
	/// @param period period of fourier function
	/// @param c coefficients of fourier function

	fourier(int num_fourier, T period, std::vector<T> c)
	: f_num_fourier(num_fourier), f_period(period), f_c(c) {}

	/// @brief copy constructer of fourier class

	fourier(const fourier<T> &copy)
	: f_num_fourier(copy.f_num_fourier), \
	f_period(copy.f_period), f_c(copy.f_c) {}

    /// @brief overloading of assignment operator for 
	/// fourier class
	fourier<T> & operator=(const fourier<T> &copy)
	{
	f_num_fourier = copy.f_num_fourier;
	f_period = copy.f_period;
	f_c = copy.f_c;
	return *this;
	}

    /// @brief update coefficients
	/// @param c coefficients of fourier function
	void update(std::vector<T> c)
	{ f_c = c;}

	/// @brief evaluate the fourier function
	/// @param t points to evaluate the fourier function
	/// @return values of the fourier function evaluated at t
	T eval(T t);

	/// @brief evaluate the fourier function
	/// @param t points to evaluate the fourier function
	/// @return values of the fourier function evaluated at t
	std::vector<T> eval(std::vector<T> t);

	/// @brief evaluate the derivative of fourier function
	/// @param t points to evaluate the derivative of fourier function
	/// @return values of the derivative of fourier function
	/// evaluated at t
	T deriv(T t);

	/// @brief evaluate the derivative of fourier function
	/// @param t points to evaluate the derivative of fourier function
	/// @return values of the derivative of fourier function
	/// evaluated at t
	std::vector<T> deriv(std::vector<T> t);

	/// @brief evaluate the nth derivative of fourier function
	/// @param n order of derivative to compute
	/// @param t points to evaluate the derivative of fourier function
	/// @return values of the nth order derivative of fourier function
	/// evaluated at t
	T nderiv(int n, T t);

	/// @brief evaluate the nth derivative of fourier function
	/// @param n order of derivative to compute
	/// @param t points to evaluate the derivative of fourier function
	/// @return values of the nth order derivative of fourier function
	/// evaluated at t
	std::vector<T> nderiv(int n, std::vector<T> t);
};
}

#include "fourier/fourier.tpp"

#endif
