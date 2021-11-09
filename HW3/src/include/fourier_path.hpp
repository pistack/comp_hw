/*!
 * @file fourier_path.hpp
 * @ingroup libpath
 * @brief headerfile for path approximated by fourier function
 * @author pistack (Junho Lee)
 * @date 2021. 11. 10.
 */

#ifndef FOURIER_PATH_H
#define FOURIER_PATH_H

#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>
#include <limits>

#include "fourier.hpp"

namespace libpath {
/// @brief Class for the path approximated by fourier function
/// \f{equation}{\label{eq:fourier_path}
/// \Psi\{\phi\{c, T\}(t)\}(t) = a + s\phi\{c, T\}(t)
/// \f}
/// ,where \f$ \phi\{c, T\}(t) \f$ is a fourier function defined in 
/// libpath::fourier ,
/// \f$ a \f$ and \f$ s \f$ are adder and scaler to match boundary
/// condition \f$ \Psi(t_0)=p_0, \Psi(t_1)=p_1 \f$, respectively.
/// @ingroup libpath
/// @param T precision should be one of float, double, long double
/// @warning If your fourier function is not vaild for 
/// the approximation of path, eval and deriv method return
/// always zero
template<typename T>
class fourier_path
{
	private:
	/// initial time
	T p_t0; 
	// final time
	T p_tf; 
	// initial coordinate of path
	T p_0;
	// final coordinate of path
	T p_f;

	/// fourier function used to approximate path
	fourier<T> p_func; 

	/// Whether or not the fourier function is valid to approximate path
	bool p_vaild; 

	/// scaler used to match initial condition
	T scale; 
	/// adder used to match initial condition
	T add; 

	/// @brief helper function for initialization
	void init_helper();

	public:

	/// @brief initalize fourier path class

	fourier_path() {}

	/// @brief initalize fourier path class
	/// @param t_init initial time
	/// @param t_fin finial time
	/// @param init initial value of path
	/// @param fin finial value of path

	fourier_path(T t_init, T t_fin,
	T init, T fin)
	: p_t0(t_init), p_tf(t_fin), p_0(init), p_f(fin)
	{}

	/// @brief initalize fourier path class
	/// @param t_init initial time
	/// @param t_fin finial time
	/// @param init initial value of path
	/// @param fin finial value of path
	/// @param fourier fourier function used to approximation of path

	fourier_path(T t_init, T t_fin,
	T init, T fin, fourier<T> fourier)
	: p_t0(t_init), p_tf(t_fin), p_0(init), p_f(fin), \
	p_func(fourier) {init_helper();}

	/// @brief copy constructor of fourier path class
	fourier_path(const fourier_path<T> &copy)
	: p_t0(copy.p_t0), p_tf(copy.p_tf), p_0(copy.p_0), p_f(copy.p_f), \
	p_func(copy.p_func), p_vaild(copy.p_vaild), scale(copy.scale), \
	add(copy.add) {}

    /// @brief overloading of assignment operator for 
	/// fourier_path class
	fourier_path<T> & operator=(const fourier_path<T> &copy)
	{
		// copy initial conditions
		p_t0 = copy.p_t0; p_tf = copy.p_tf; p_0 = copy.p_0; p_f = copy.p_f; 
		// copy fourier function, scaler and adder
		p_func = copy.p_func; p_vaild = copy.p_vaild; 
		scale = copy.scale; add = copy.add;
		return *this;
	}

	/// @brief update fourier function
	/// It also updates the validity of fourier function
	/// @param fourier fourier function to update
	void update(fourier<T> fourier)
	{
		p_func = fourier;
		init_helper();
	}

    /// @brief check whether or not the fourier function is valid to
	/// approximate path
	/// @return vaildity of fourier function
	bool is_vaild()
	{
		return p_vaild;
	}

	/// @brief get initial and final time of path
	/// @return tuple of initial and final time of path
	std::tuple<T, T> get_endtimes()
	{
		return std::make_tuple(p_t0, p_tf);
	}

	/// @brief get adder
	/// @return adder
	T get_adder()
	{
		return add;
	}

	/// @brief get scaler
	/// @return scaler
	T get_scaler()
	{
		return scale;
	}

	/// @brief evaluate the path
	/// @param t points to evaluate the path
	/// @return the path evaluated at t
	T eval(T t)
	{
		return scale*p_func.eval(t)+add;
	}

	/// @brief evaluate the path
	/// @param t points to evaluate the path
	/// @return the path evaluated at t
	std::vector<T> eval(std::vector<T> t);

	/// @brief evaluate the derivative of path
	/// @param t points to evaluate
	/// @return the derivative of path evaluated at t
	T deriv(T t)
	{
		return scale*p_func.deriv(t);
	}

	/// @brief evaluate the path
	/// @param t points to evaluate the path
	/// @return the path evaluated at t
	std::vector<T> deriv(std::vector<T> t);

	/// @brief evaluate the nth derivative of path
	/// @param n order of derivative to compute
	/// @param t points to evaluate
	/// @return the nth order derivative of path evaluated at t
	T nderiv(int n, T t)
	{
		return scale*p_func.nderiv(n, t);
	}

	/// @brief evaluate the nth derivative of path
	/// @param n order of derivative to compute
	/// @param t points to evaluate
	/// @return the nth order derivative of path evaluated at t
	std::vector<T> nderiv(int n, std::vector<T> t);
};
}

#include "fourier_path/fourier_path.tpp"

#endif