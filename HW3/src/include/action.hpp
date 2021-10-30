/*!
 * @ingroup libfourier
 * @file action.hpp
 * @brief header file for evaluation of the action
 * @author pistack (Junho Lee)
 * @date 2021. 10. 31.
 */

#ifndef ACTION_H
#define ACTION_H

#include <cmath>
#include <limits>
#include <vector>
#include "fourier_path.hpp"

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
/// @see [J. ACM 16, 3 (July 1969), 
/// 483–495.](https://dl.acm.org/doi/10.1145/321526.321537)
/// @ingroup libfourier
template<typename T, typename Lag>
class action
{
	private:
	T atol; // absoulte tolerance
	
	std::vector<fourier_path<T>> path_action; // path

	bool vaildity=false; // vaildity of path

	T D_tol; // tolerance of D factor

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
	/// @param D previous D value 
	T eval_helper(T left, T mid, T right,
	T fleft, T fmid, T fright, 
	T integral, T D);

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
	/// by adpative simpson's method
	/// @return action of given path

	T eval();
};

#include "action/action.tpp"

#endif