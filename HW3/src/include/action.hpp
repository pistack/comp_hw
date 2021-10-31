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
/// @see 
/// @ingroup libfourier
template<typename T, typename Lag>
class action
{
	private:
	// node and weight for G15-K31 Gauss-Kronrod Quadrature
	const std::vector<T> nodes{
		9.980022986933970602851728401522712e-01, 9.879925180204854284895657185866126e-01,
		9.677390756791391342573479787843372e-01, 9.372733924007059043077589477102095e-01,
		8.972645323440819008825096564544959e-01, 8.482065834104272162006483207742169e-01,
		7.904185014424659329676492948179473e-01, 7.244177313601700474161860546139380e-01,
		6.509967412974169705337358953132747e-01, 5.709721726085388475372267372539106e-01,
		4.850818636402396806936557402323506e-01, 3.941513470775633698972073709810455e-01,
		2.991800071531688121667800242663890e-01, 2.011940939974345223006283033945962e-01,
		1.011420669187174990270742314473923e-01	
	};
	const std::vector<T> weight_gauss{
		3.075324199611726835462839357720442e-02, 7.036604748810812470926741645066734e-02,
		1.071592204671719350118695466858693e-01, 1.395706779261543144478047945110283e-01,
		1.662692058169939335532008604812088e-01, 1.861610000155622110268005618664228e-01,
		1.984314853271115764561183264438393e-01, 2.025782419255612728806201999675193e-01
	};
  const std::vector<T> weight_kronrod{
	  5.377479872923348987792051430127650e-03, 1.500794732931612253837476307580727e-02,
	  2.546084732671532018687400101965336e-02, 3.534636079137584622203794847836005e-02,
	  4.458975132476487660822729937327969e-02, 5.348152469092808726534314723943030e-02,
	  6.200956780067064028513923096080293e-02, 6.985412131872825870952007709914748e-02,
	  7.684968075772037889443277748265901e-02, 8.308050282313302103828924728610379e-02,
	  8.856444305621177064727544369377430e-02, 9.312659817082532122548687274734572e-02,
	  9.664272698362367850517990762758934e-02, 9.917359872179195933239317348460313e-02,
	  1.007698455238755950449466626175697e-01, 1.013300070147915490173747927674925e-01
	};
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

	/// @brief evaluate lagranian at given t
	/// @param t time at which evaluate lagrangian
	/// @return value of the lagranian at given t
	std::vector<T> eval_lagranian(std::vector<T> t);

	/// @brief helper function for action evaluation
	/// (G15, K31) Gauss–Kronrod quadrature method
    /// @param left left end point of interval
	/// @param right right end point of interval
	/// @param D previous |G15-K31| value without scaling
	T eval_helper(T left, T right, T D);

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
	/// by (G7, K15) Gauss–Kronrod quadrature method
	/// @return action of given path

	T eval();
};

#include "action/action.tpp"

#endif