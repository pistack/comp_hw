/*!
 * @file hw3.hpp
 * @brief headerfile for homework3 of Computer1 class in Yonsei University
 * Minimize the action by Markov Chain Monte Carlo Method 
 * to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 10.
 */

#ifndef HW3_H
#define HW3_H
#endif

#include <cmath>
#include <random>
#include <tuple>
#include <vector>

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
 * @brief evaluates sum of sine function with 
 * period 4*(t_max-t_min) weighted by coefficents and its derivative
 * at given time t.
 * @param t time
 * @param c coefficiets
 * @param num_sine number of sine function to add
 * @return tuple of value and derivatives of the sum of sine function
 */

std::tuple<std::vector<double>, std::vector<double>>
sum_of_sine(std::vector<double> &t, std::vector<double> c, int num_sine);

/*!
 * @brief randomly move initial guess by step
 * @param init_guess initial guess
 * @param step step size\
 * @param gen random number generator (assume mt19937)
 * @param dist distribution
 * @return moved initial guess by step
 */

std::vector<double>
move_step(std::vector<double> &init_guess, double step,
	  std::mt19937 &gen, std::uniform_real_distribution<double> &dist);

/*!
 * @brief evaluates the action of given path
 * @param t time
 * @param zeta zeta part of path
 * @param deriv_zeta derivative of path, zeta part
 * @param theta theta part of path
 * @param deriv_theta derivative of path, theta part
 * @return the action of given path
 */

double eval_action(std::vector<double> &t, std::vector<double> zeta,
		   std::vector<double> deriv_zeta, std::vector<double> theta,
		   std::vector<double> deriv_theta);

/*!
 * @brief HW3: Solve Kepler problem via Markov Chain Monte Carlo Method
 *  see HW3.pdf for further detail.
 * @param zeta_min minimum value of zeta, 
    for constraint motion 0.5 < zeta_min < 1
 * @param t0 initial time
 * @param n number of points to evaluate
 * @param num_sine number of sine functions to guess
 * @param num_iter number of iteration
 * @param step step size
 * @param lambda paramter to adapt step size
 * @param gen random number generator (assume mt19937)
 * @param dist distribution
 * @return tuple of number of actural moves, minimum action, 
 *  time and path(zeta, theta)
 */

std::tuple<int, double,
	   std::vector<double>, std::vector<double>, std::vector<double>>
HW3(double zeta_min, double t0,
    int n, int num_sine, int num_iter, double step, double lambda,
    std::mt19937 &gen, std::uniform_real_distribution<double> &dist);
