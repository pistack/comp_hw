/*!
 * @file fourier.hpp
 * @brief headerfile for fourier function
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#ifndef FOURIER_H
#define FOURIER_H
#endif

#include <vector>

/*!
 * @brief Class for sum of the sine and cosine function
 * weighted by coefficients
 */

class fourier
{
	private:
        const double pi = 3.141592653589793; // define pi 
	int f_num_fourier; ///< number of sine and cosine function to add
	double f_period; ///< the period of fourier function 
	std::vector<double> f_c; ///< coefficients

	public:
	/*!
	 * @brief initize fourier class
	 * @param num_fourier number of sine and cosine function to add
	 * @param period period of fourier function
	 * @param c coeffcients of fourier function
	 */

	void init(int num_fourier, double period, std::vector<double> &c);

	/*!
	 * @brief evaluate the fourier function
	 * @param t points to evaluate the fourier function
	 * @return values of the fourier function evaluated at t
	 */

	std::vector<double> eval(std::vector<double> &t);

	/*!
	 * @brief evaluate the derivative of fourier function
	 * @param t points to evaluate the derivative of fourier function
	 * @return values of the derivative of fourier function
	 * evaluated at t
	 */

	std::vector<double> deriv(std::vector<double> &t);
};

/*!
 * @brief Class for the path approximated by fourier function
 */

class fourier_path
{
	private:
	///< fourier function used to approximate path
	fourier p_func; 
	double p_init; ///< initial value of path
	double p_final; ///< final value of path
	double scale; ///< scaler used to match initial condtion
	double add; ///< adder used to match initial condition

	public:

	/*!
	 * @brief initize fourier path class
	 * @param t_init initial time
	 * @param t_fin finial time
	 * @param init initial value of path
	 * @param fin finial value of path
	 * @param period period of fourier function
	 * @param c coefficients of fourier function
	 */

	void init(double t_init, double t_fin,
	double init, double fin, double period,
	std::vector<double> &c);

	/*!
	 * @brief evaluate the path
	 * @param t points to evaluate the path
	 * @return the path evaluated at t
	 */

	std::vector<double> eval(std::vector<double> &t);

	/*!
	 * @brief evaluate the derivative of path
	 * @param t points to evaluate
	 * @return the derivative of path evaluated at t
	 */

	std::vector<double> deriv(std::vector<double> &t);
};


