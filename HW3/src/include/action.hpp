/*!
 * @ingroup libpath
 * @file action.hpp
 * @brief header file for evaluation of the action
 * @author pistack (Junho Lee)
 * @date 2021. 11. 10.
 */

#ifndef ACTION_H
#define ACTION_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "node_weight_table.hpp"

namespace libpath{

template<typename T>
constexpr T h_pi = T(std::acos(0)); ///< half pi

/// @brief class which computes action
/// @param T precision should be one of
/// float, double and long double
/// @param Path type of path
/// should be one of fourier_path or bezier_path
/// @param Lag lagranian of action
/// functor class which has
/// time, path and derivative of path
/// as variable and it returns
/// value of lagranian at given time
/// @warning If you give invaild path, then
/// eval method will return always zero.
///
/// For gauss-kronrod quadrature method,
/// @see [Math. Comp. 22 (1968), 847-856.]
/// (https://doi.org/10.1090/S0025-5718-68-99866-9) 
/// @see [Commun. ACM 16, 11 (Nov. 1973), 694–699.]
/// (https://doi-org-ssl.access.yonsei.ac.kr/10.1145/355611.362543)
/// @see [ACM Comput. Surv. 44, 4, Article 22 (August 2012)]
/// (https://doi.org/10.1145/2333112.2333117)
///
/// For tanh-sinh quadrature method, 
/// @see [Publ. RIMS, Kyoto Univ. 9 (1974), 721-741]
/// (https://doi.org/10.2977%2Fprims%2F1195192451)
/// @see [Publ. RIMS, Kyoto Univ. 41 (2005), 897-935]
/// (https://doi.org/10.2977%2Fprims%2F1145474600)
/// @see [David H. Bailey, Tanh-Sinh High-Precision Quadrature]
/// (https://www.davidhbailey.com/dhbpapers/dhb-tanh-sinh.pdf)
/// @ingroup libpath
template<typename T, typename Path, typename Lag>
class action
{
	private:

	T atol; // absoulte tolerance
	
	std::vector<Path> path_action; // path

	bool vaildity=false; // vaildity of path

	const Lag lag; // lagrangian

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
	/// @param[out] integral integration value
	/// @param[out] e estimated error
	template<typename Gau_Kron>
	void eval_helper(T left, T right, T D, T D_tol, T &integral, T &e);

	/// @brief evaluate the action of given path
	/// by \f$ (G_{(n-1)/2}, K_n) \f$ Gauss–Kronrod quadrature method
	/// @param left left end points of interval
	/// @param right right end points of interval
	/// @param n order of gauss-kronrod quadrature
	/// currently only supports N=15, 21, 31, 41, 51, 61
	/// @param[out] e estimated error
	/// @return action of given path

	T eval_quadgk(T left, T right, int n, T &e);

	/// @brief evaluate the action of given path
	/// tanh-sinh quadrature method
	/// @param left left end points of interval
	/// @param right right end points of interval
	/// @param max_depth maximum depth 
	/// step size \f$ h \f$ would be 
	/// \f$ h \geq 2^{-\mathrm{depth}_{max}} \f$.
	/// @param[out] e estimated error 
	/// @return action of given path

	T eval_qthsh(T left, T right, int max_depth, T &e);

	public:
    /// @brief initialize action class
	action(){}

	/// @brief initialize action class
	/// @param path path

	action(std::vector<Path> path)
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
	std::vector<Path> path)
	: atol(abs_tol), path_action(path) 
	{check_vaild();}

	/// @brief copy constructer of action class
	action(const action<T, Path, Lag> &copy)
	: atol(copy.atol), path_action(copy.path_action), \
	vaildity(copy.vaildity)
	{}

    /// @brief overloading of assignment operator for 
	/// action class
	action<T, Path, Lag> & operator=(const action<T, Path, Lag> &copy)
	{
		 atol = copy.atol;
		 path_action = copy.path_action;
		 vaildity = copy.vaildity;
		 return *this;
	}

	/// @brief update path
	/// @param path path to update
	void update(std::vector<Path> path)
	{path_action = path; check_vaild();}

    /// @brief update absolute tolerance
	/// @param abs_tol absolute tolerance
	void update(T abs_tol)
	{atol=abs_tol;}

	/// @brief update absolute tolerance and path
	/// @param path path to update
	/// @param abs_tol absolute tolerance
	void update(std::vector<Path> path,
	T abs_tol)
	{atol=abs_tol; path_action = path; check_vaild();}

	/// @brief check vaildity of path
	/// @return vaildity of path
	bool is_vaild()
	{return vaildity;}

	/// @brief evaluate the action of given path
	/// by default method:
	/// (G15, K31) Gauss–Kronrod quadrature method
	/// @param[out] e estimated error 
	/// @return action of given path

	T eval(T &e);

	/// @brief evaluate the action of given path
	/// @param method numerical integration method
	/// - method 0: Gauss-Kronrod quadrature method
	/// - method 1: Tanh-Sinh quadrature method
	/// @param n  
	/// - order of Gauss-Kronrod quadrature if method equals to 0,
	///   currently, only support n=15, 21, 31, 41, 51, 61.
	/// - maximum depth if method equals to 1.
	/// @param[out] e estimated error
	/// @return action of given path

	T eval(int method, int n, T &e);
};
}

#include "action/action.tpp"

#endif