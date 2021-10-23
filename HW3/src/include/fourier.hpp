/*!
 * @file fourier.hpp
 * @ingroup libfourier
 * @brief headerfile for fourier function
 * @author pistack (Junho Lee)
 * @date 2021. 10. 18.
 */

#ifndef FOURIER_H
#define FOURIER_H

#include <vector>

/// @brief Class for sum of the sine and cosine function
/// @ingroup libfourier
/// weighted by coefficients
class fourier
{
	private:
        const double pi = 3.141592653589793; // pi 
		int f_num_fourier; // number of sine and cosine function to add
		double f_period; // the period of fourier function 
		std::vector<double> f_c; // coefficients

	public:


	/// @brief initialize fourier class
	fourier() {}

	/// @brief initialize fourier class
	/// @param num_fourier number of sine and consine function to add
	/// @param period period of fourier function
	/// @param c coefficients of fourier function

	fourier(int num_fourier, double period, std::vector<double> c)
	: f_num_fourier(num_fourier), f_period(period), f_c(c) {}

	/// @brief copy constructer of fourier class

	fourier(const fourier &copy)
	: f_num_fourier(copy.f_num_fourier), \
	f_period(copy.f_period), f_c(copy.f_c) {}

	fourier& operator=(const fourier &copy);

    /// @brief update coefficients
	/// @param c coefficients of fourier function
	void update(std::vector<double> c);

	/// @brief evaluate the fourier function
	/// @param t points to evaluate the fourier function
	/// @return values of the fourier function evaluated at t
	double eval(double t);

	/// @brief evaluate the fourier function
	/// @param t points to evaluate the fourier function
	/// @return values of the fourier function evaluated at t
	std::vector<double> eval(std::vector<double> t);


	/// @brief evaluate the derivative of fourier function
	/// @param t points to evaluate the derivative of fourier function
	/// @return values of the derivative of fourier function
	/// evaluated at t
	double deriv(double t);

	/// @brief evaluate the derivative of fourier function
	/// @param t points to evaluate the derivative of fourier function
	/// @return values of the derivative of fourier function
	/// evaluated at t
	std::vector<double> deriv(std::vector<double> t);

	/// @brief evaluate the nth derivative of fourier function
	/// @param n order of derivative to compute
	/// @param t points to evaluate the derivative of fourier function
	/// @return values of the nth order derivative of fourier function
	/// evaluated at t
	double nderiv(int n, double t);

	/// @brief evaluate the nth derivative of fourier function
	/// @param n order of derivative to compute
	/// @param t points to evaluate the derivative of fourier function
	/// @return values of the nth order derivative of fourier function
	/// evaluated at t
	std::vector<double> nderiv(int n, std::vector<double> t);
};

#endif


