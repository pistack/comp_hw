/*!
 * @file hw3.hpp
 * @ingroup hw3
 * @brief headerfile for homework3 of Computer1 class in Yonsei University
 * Minimize the action by Markov Chain Monte Carlo Method 
 * to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 10.
 */

#ifndef HW3_H
#define HW3_H

#include <tuple>
#include <vector>
#include <random>

const double pi = 3.141592653589793; ///<define pi

/*!
 * @brief randomly move initial guess by at most step
 * @param init_guess initial guess
 * @param step step size
 * @param gen random number generator (assume mt19937)
 * @param dist distribution
 * @return moved initial guess by at most step
 * @ingroup hw3
 */

std::vector<double>
move_step(std::vector<double> init_guess, double step,
	  std::mt19937 &gen, std::uniform_real_distribution<double> &dist);

/*!
 * @brief HW3: Solve Kepler problem via Markov Chain Monte Carlo Method
 * @param t0 initial time
 * @param zeta_min minimum value of zeta, 
    for constraint motion 0.5 < zeta_min < 1
 * @param atol absolute tolerance of action
 * @param rtol relative tolerance of action
 * @param num_fourier number of sine and cosine functions to guess
 * @param num_eval number of points to evaluate minimum path
 * @param max_iter maximum number of iteration
 * @param conv_atol convergence criterion: changes of abs value of minimum action
 * @param conv_rtol convergence criterion: changes of rel value of minimum action
 * @param step step size
 * @param lambda paramter to adapt step size
 * @param gen random number generator (assume mt19937)
 * @param dist distribution
 * @return tuple of the number of actural moves, minimum action, 
 *  time and path(zeta, theta) evaluated at time
 * @ingroup hw3
 */

std::tuple<int, double, 
std::vector<double>, std::vector<double>, std::vector<double>>
HW3(double t0, double zeta_min, double atol, double rtol, 
int num_fourier, int num_eval,
int max_iter, double conv_atol, double conv_rtol, 
double step, double lambda,
std::mt19937 &gen, std::uniform_real_distribution<double> &dist);

#endif
