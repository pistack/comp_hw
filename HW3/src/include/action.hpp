/*!
 * @ingroup libfourier
 * @file action.hpp
 * @brief header file for evaluation of the action
 * @author pistack (Junho Lee)
 * @date 2021. 10. 24.
 */

#ifndef ACTION_H
#define ACTION_H

#include <vector>
#include "fourier_path.hpp"

/// @brief class to compute action
/// @warning If you give invaild path, then
/// eval method will return always zero
/// @ingroup libfourier
class action
{
	private:
	const double MAXDEPTH = 30; // maximum depth of recurrsion
	double atol, rtol; // abs and relative tol
	double (*lagranian)(double, std::vector<double>,
	std::vector<double>);
	std::vector<fourier_path> path_action; // path
	// lagranian of action
	bool vaildity=false; // vaildity of path

	/// @brief checks the vaildity of path
	void check_vaild();

	/// @brief evaluate lagranian at given t
	/// @param t time to evaluate lagrangian
	/// @return value of the lagranian at given t
	double eval_lagranian(double t);

    /// @brief helper function for action evaluation
    /// @param left left end point of interval
	/// @param mid mid point of interval
	/// @param right right end point of interval
	/// @param fleft value of lagrangian at left end point
	/// @param fmid value of lagrangian at mid pooint
	/// @param fright value of lagrangian at right end point
	/// @param integral integrated value
	/// @param tol tolerance
	/// @param depth recurrsion depth
	double eval_helper(double left, double mid, double right,
	double fleft, double fmid, double fright, 
	double integral, double tol, int depth);

	public:
    /// @brief initialize action class
	action(){}

	/// @brief initialize action class
	/// @param abs_tol absoulte tolerance
	/// @param rel_tol relative tolerance
	/// @param lag lagranian of action
	/// - first param:  t
	/// - type of first param: double 
	/// - second param: double path at t
	/// - type of second param: vector<double>
	/// - third param: derivative of path at t
	/// - type of third param: vector<double>
	/// - type of return value: double

	action(double abs_tol, double rel_tol, 
	double (*lag)(double, std::vector<double>, 
	std::vector<double>))
	: atol(abs_tol), rtol(rel_tol), \
	lagranian(lag)
	{}

	/// @brief initialize action class
	/// @param abs_tol absoulte tolerance
	/// @param rel_tol relative tolerance
	/// @param path path
	/// @param lag lagranian of action
	/// - first param:  t
	/// - type of first param: double 
	/// - second param: double path at t
	/// - type of second param: vector<double>
	/// - third param: derivative of path at t
	/// - type of third param: vector<double>
	/// - type of return value: double

	action(double abs_tol, double rel_tol,
	std::vector<fourier_path> path, 
	double (*lag)(double, std::vector<double>, 
	std::vector<double>))
	: atol(abs_tol), rtol(rel_tol), \
	lagranian(lag), path_action(path) 
	{check_vaild();}

	/// @brief copy constructer of action class
	action(const action &copy)
	: atol(copy.atol), rtol(copy.rtol), \
	lagranian(copy.lagranian), path_action(copy.path_action), \
	vaildity(copy.vaildity)
	{}

	action & operator=(const action &copy);

	/// @brief update path
	/// @param path path to update
	void update(std::vector<fourier_path> path);

	/// @brief check vaildity of path
	/// @return vaildity of path
	bool is_vaild();

	/// @brief evaluate the action of given path
	/// by adpative simpson's method
	/// @return action of given path

	double eval();
};

#endif