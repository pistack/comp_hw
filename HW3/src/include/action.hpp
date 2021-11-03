/*!
 * @ingroup libfourier
 * @file action.hpp
 * @brief header file for evaluation of the action
 * @author pistack (Junho Lee)
 * @date 2021. 11. 3.
 */

#ifndef ACTION_H
#define ACTION_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "fourier_path.hpp"
#include "node_weight_table.hpp"

namespace libfourier{
/// @brief class which computes action
/// @param T precision should be one of
/// float, double and long double
/// @param Lag lagranian of action
/// functor class which has
/// time, path and derivative of path
/// as variable and it returns
/// value of lagranian at given time
/// @warning finial time of path should be greater than
/// initial time of it.
/// @warning If you give invaild path, then
/// eval method will return always zero.
/// @see for gauss-kronrod quadrature method  
/// @see for tanh-sinh method
/// @ingroup libfourier
template<typename T, typename Lag>
class action
{
	private:

	T atol; // absoulte tolerance
	
	std::vector<fourier_path<T>> path_action; // path

	bool vaildity=false; // vaildity of path

	/// @brief checks the vaildity of path
	void check_vaild();

	/// @brief evaluate lagrangian at given t
	/// @param t time at which evaluate lagrangian
	/// @return value of the lagrangian at given t
	T eval_lagrangian(T t);

	/// @brief evaluate lagrangian at given t
	/// @param t time at which evaluate lagrangian
	/// @return value of the lagrangian at given t
	std::vector<T> eval_lagrangian(std::vector<T> t);

	/// @brief helper function for action evaluation by
	/// \f$ (G_{(n-1)/2}, Kn) \f$ Gauss–Kronrod quadrature method
	/// @param Gau_Kron table for gauss_kronrod node and weights
    /// @param left left end point of interval
	/// @param right right end point of interval
	/// @param D previous \f$ |G_{(n-1)/2} - K_n| \f$. value without scaling
	/// @param D_tol tolerance for D
	template<typename Gau_Kron>
	T eval_helper(T left, T right, T D, T D_tol);

	/// @brief evaluate the action of given path
	/// by \f$ (G_{(n-1)/2}, K_n) \f$ Gauss–Kronrod quadrature method
	/// @param left left end points of interval
	/// @param right right end points of interval
	/// @param n order of gauss-kronrod quadrature
	/// currently only supports N=15, 21, 31, 41, 51, 61
	/// @return action of given path

	T eval_quadgk(T left, T right, int n);

	/// @brief evaluate the action of given path
	/// tanh-sinh quadrature method
	/// @param left left end points of interval
	/// @param right right end points of interval
	/// @param max_depth maximum depth 
	/// step size \f$ h \f$ would be 
	/// \f$ h \geq 2^{-\mathrm{depth}_{max}} \f$. 
	/// @return action of given path

	T eval_qthsh(T left, T right, int max_depth);

	public:
    /// @brief initialize action class
	action(){}

	/// @brief initialize action class
	/// @param path path

	action(std::vector<fourier_path<T>> path)
	: path_action(path)
	{check_vaild();}

	/// @brief initialize action class
	/// @param abs_tol absoulte tolerance

	action(T abs_tol)
	: atol(abs_tol)
	{}

	/// @brief initialize action class
	/// @param abs_tol absoulte tolerance
	/// @param path path

	action(T abs_tol,
	std::vector<fourier_path<T>> path)
	: atol(abs_tol), path_action(path) 
	{check_vaild();}

	/// @brief copy constructer of action class
	action(const action<T, Lag> &copy)
	: atol(copy.atol), path_action(copy.path_action), \
	vaildity(copy.vaildity)
	{}

	action<T, Lag> & operator=(const action<T, Lag> &copy);

	/// @brief update path
	/// @param path path to update
	void update(std::vector<fourier_path<T>> path);

    /// @brief update absolute tolerance
	/// @param abs_tol absolute tolerance
	void update(T abs_tol);

	/// @brief update absolute tolerance and path
	/// @param path path to update
	/// @param abs_tol absolute tolerance
	void update(std::vector<fourier_path<T>> path,
	T abs_tol);

	/// @brief check vaildity of path
	/// @return vaildity of path
	bool is_vaild();

	/// @brief evaluate the action of given path
	/// by default method:
	/// (G15, K31) Gauss–Kronrod quadrature method
	/// @return action of given path

	T eval();

	/// @brief evaluate the action of given path
	/// @param method numerical integration method
	/// - method 0: Gauss-Kronrod quadrature method
	/// - method 1: Tanh-Sinh quadrature method
	/// @param n  
	/// - order of Gauss-Kronrod quadrature if method equals to 0,
	/// - maximum depth if method equals to 1.
	/// @return action of given path

	T eval(int method, int n);
};
}

#include "action/action.tpp"

#endif