/*!
 * @file kepler.hpp
 * @brief headerfile for evaluate the kepler action
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#ifndef KEPLER_H
#define KEPLER_H
#endif

#include <vector>

/*!
 * @brief evaluates the action of given path
 * @param t time
 * @param zeta_init initial value of zeta
 * @param zeta_finial final value of zeta
 * @param theta_init initial value of theta
 * @param theta_finial final value of theta
 * @param c_zeta fourier coefficients for zeta
 * @param c_theta fourier coefficients for theta
 * @return the action of given path
 */

double eval_action(std::vector<double> &t,
		   double zeta_init, double zeta_finial,
		   double theta_init, double theta_finial,
		   std::vector<double> c_zeta,
		   std::vector<double> c_theta);
