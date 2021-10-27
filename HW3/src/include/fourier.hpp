/*!
 * @file fourier.hpp
 * @ingroup libfourier
 * @brief headerfile for fourier function
 * @author pistack (Junho Lee)
 * @date 2021. 10. 28.
 */

#ifndef FOURIER_H
#define FOURIER_H

#include <algorithm>
#include <cmath>
#include <vector>

/// @brief Class for sum of the sine and cosine function
/// weighted by coefficients
/// class is undefined when type T is not one of
/// float, double, long double
/// @ingroup libfourier

template<typename T>
class fourier
{
	private:
        const T pi = \
		3.1415926535897932384626433832795028841971693993751058209749445923078164062; // pi 
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

	fourier<T> & operator=(const fourier<T> &copy);

    /// @brief update coefficients
	/// @param c coefficients of fourier function
	void update(std::vector<T> c);

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

#include "fourier/fourier.tpp"

#endif