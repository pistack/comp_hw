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

/// @brief Class for sum of the sine and cosine function
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
	/// @param num_fourier number of sine and consine function to add
	/// @param period period of fourier function
	/// @param c coefficients of fourier function

	fourier() {}

	fourier(int num_fourier, double period, std::vector<double> c)
	: f_num_fourier(num_fourier), f_period(period), f_c(c) {}

	/// @brief copy constructer for fourier class

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
	std::vector<double> eval(std::vector<double> t);


	/// @brief evaluate the derivative of fourier function
	/// @param t points to evaluate the derivative of fourier function
	/// @return values of the derivative of fourier function
	/// evaluated at t
	double deriv(double t);
	std::vector<double> deriv(std::vector<double> t);

	/// @brief evaluate the nth derivative of fourier function
	/// @param n order of derivative to compute
	/// @param t points to evaluate the derivative of fourier function
	/// @return values of the nth order derivative of fourier function
	/// evaluated at t
	double nderiv(int n, double t);
	std::vector<double> nderiv(int n, std::vector<double> t);
};

/// @brief Class for the path approximated by fourier function
/// @warning If your fourier function is not vaild for 
/// the approximation of path, eval and deriv method return
/// always zero
class fourier_path
{
	private:
	/// initial time
	double p_t0; 
	// final time
	double p_tf; 
	// initial coordinate of path
	double p_0;
	// final coordinate of path
	double p_f;

	/// fourier function used to approximate path
	fourier p_func; 

	/// Whether or not the fourier function is valid to approximate path
	bool p_vaild; 

	/// scaler used to math initial condition
	double scale; 
	/// adder used to match initial condition
	double add; 

	/// @brief helper function for initialization
	void init_helper();

	public:

	/// @brief initalize fourier path class
	/// @param t_init initial time
	/// @param t_fin finial time
	/// @param init initial value of path
	/// @param fin finial value of path
	/// @param fourier fourier function used to approximation of path

	fourier_path() {}

	fourier_path(double t_init, double t_fin,
	double init, double fin, fourier fourier)
	: p_t0(t_init), p_tf(t_fin), p_0(init), p_f(fin), \
	p_func(fourier) {init_helper();}

	/// @brief copy constructor of fourier path class
	fourier_path(const fourier_path &copy)
	: p_t0(copy.p_t0), p_tf(copy.p_tf), p_0(copy.p_0), p_f(copy.p_f), \
	p_func(copy.p_func), p_vaild(copy.p_vaild), scale(copy.scale), \
	add(copy.add) {}

	fourier_path& operator=(const fourier_path &copy);

	/// @brief update fourier function
	/// It also updates the validity of fourier function
	/// @param fourier fourier function to update
	void update(fourier fourier);

    /// @brief check whether or not the fourier function is valid to
	/// approximate path
	/// @return vaildity of fourier function
	bool is_vaild();

	/// @brief evaluate the path
	/// @param t points to evaluate the path
	/// @return the path evaluated at t
	double eval(double t);
	std::vector<double> eval(std::vector<double> t);

	/// @brief evaluate the derivative of path
	/// @param t points to evaluate
	/// @return the derivative of path evaluated at t
	double deriv(double t);
	std::vector<double> deriv(std::vector<double> t);

	/// @brief evaluate the nth derivative of path
	/// @param n order of derivative to compute
	/// @param t points to evaluate
	/// @return the nth order derivative of path evaluated at t
	double nderiv(int n, double t);
	std::vector<double> nderiv(int n, std::vector<double> t);
};


