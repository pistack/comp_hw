/*!
 * @ingroup libfourier
 * @file action.hpp
 * @brief header file for evaluation of the action
 * @author pistack (Junho Lee)
 * @date 2021. 10. 28.
 */

#ifndef ACTION_H
#define ACTION_H

#include <cmath>
#include <cerrno>
#include <vector>
#include "fourier_path.hpp"

/// @brief class which computes action
/// Only defined when type T is one of float, double
/// and long double.
/// @warning If you give invaild path, then
/// eval method will return always zero
/// @ingroup libfourier
template<typename T>
class action
{
	private:
	const int MAXDEPTH = 30; // maximum depth of recurrsion
	T atol, rtol; // absoulte and relative tol

	// lagranian of action
	T (*lagranian)(T, std::vector<T>,
	std::vector<T>);

	std::vector<fourier_path<T>> path_action; // path

	bool vaildity=false; // vaildity of path

	/// @brief checks the vaildity of path
	void check_vaild();

	/// @brief evaluate lagranian at given t
	/// @param t time at which evaluate lagrangian
	/// @return value of the lagranian at given t
	T eval_lagranian(T t);

    /// @brief helper function for action evaluation
    /// @param left left end point of interval
	/// @param mid mid point of interval
	/// @param right right end point of interval
	/// @param fleft value of lagrangian at left end point
	/// @param fmid value of lagrangian at mid point
	/// @param fright value of lagrangian at right end point
	/// @param integral integrated value
	/// @param tol tolerance
	/// @param depth recurrsion depth
	T eval_helper(T left, T mid, T right,
	T fleft, T fmid, T fright, 
	T integral, T tol, int depth);

	public:
    /// @brief initialize action class
	action(){}

	/// @brief initialize action class
	/// @param abs_tol absoulte tolerance
	/// @param rel_tol relative tolerance
	/// @param lag lagranian of action
	/// - first param:  t
	/// - type of first param: T 
	/// - second param: T path at t
	/// - type of second param: vector<T>
	/// - third param: derivative of path at t
	/// - type of third param: vector<T>
	/// - type of return value: T

	action(T abs_tol, T rel_tol, 
	T (*lag)(T, std::vector<T>, 
	std::vector<T>))
	: atol(abs_tol), rtol(rel_tol), \
	lagranian(lag)
	{}

	/// @brief initialize action class
	/// @param abs_tol absoulte tolerance
	/// @param rel_tol relative tolerance
	/// @param path path
	/// @param lag lagranian of action
	/// - first param:  t
	/// - type of first param: T 
	/// - second param: T path at t
	/// - type of second param: vector<T>
	/// - third param: derivative of path at t
	/// - type of third param: vector<T>
	/// - type of return value: T

	action(T abs_tol, T rel_tol,
	std::vector<fourier_path<T>> path, 
	T (*lag)(T, std::vector<T>, 
	std::vector<T>))
	: atol(abs_tol), rtol(rel_tol), \
	lagranian(lag), path_action(path) 
	{check_vaild();}

	/// @brief copy constructer of action class
	action(const action<T> &copy)
	: atol(copy.atol), rtol(copy.rtol), \
	lagranian(copy.lagranian), path_action(copy.path_action), \
	vaildity(copy.vaildity)
	{}

	action<T> & operator=(const action<T> &copy);

	/// @brief update path
	/// @param path path to update
	void update(std::vector<fourier_path<T>> path);

	/// @brief check vaildity of path
	/// @return vaildity of path
	bool is_vaild();

	/// @brief evaluate the action of given path
	/// by adpative simpson's method
	/// @return action of given path

	T eval();
};

#include "action/action.tpp"

#endif