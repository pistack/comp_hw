/*!
 * @file basic.hpp
 * @brief headerfile for basic function
 * to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#ifndef BASIC_H
#define BASIC_H
#endif

#include <tuple>
#include <vector>
#include <random>

const double pi = 3.141592653589793; // define pi 

/*!
 * @brief scale vector by scaler and add vector by adder
 * @param v vector to scale and add
 * @param scale scaler
 * @param add adder
 * @return scale and added vector (scale*v+add)
 */

std::vector<double>
scale_and_add_vector(std::vector<double> &v, double scale, double add);

/*!
 * @brief evaluates sum of sine and cosine function with 
 * period 2*(t_max-t_min) weighted by coefficents and its derivative
 * at given time t.
 * @param t time
 * @param c coefficiets
 * @param num_fourier number of sine and cosine function to add
 * @return tuple of value and derivatives of 
 * the sum of sine and cosine function
 */

std::tuple<std::vector<double>, std::vector<double>>
sum_of_fourier(std::vector<double> &t, std::vector<double> c, int num_fourier);

/*!
 * @brief randomly move initial guess by at most step
 * @param init_guess initial guess
 * @param step step size
 * @param gen random number generator (assume mt19937)
 * @param dist distribution
 * @return moved initial guess by step
 */

std::vector<double>
move_step(std::vector<double> &init_guess, double step,
	  std::mt19937 &gen, std::uniform_real_distribution<double> &dist);

/*!
 * @brief evaluates path and its derivative
 * @param t time
 * @param init initial value of the path
 * @param fin final value of the path
 * @param c fourier coefficients for path
 * @return the tuple of path and its derivative
 */

std::tuple<std::vector<double>, std::vector<double>>
eval_path(std::vector<double> &t, double init, double finial,
	  std::vector<double> &c);
