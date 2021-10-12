/*!
 * @file kepler.hpp
 * @brief header file for evaluation of the kepler action
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#ifndef KEPLER_H
#define KEPLER_H
#endif

#include <vector>

/*!
 * @brief evaluates the action of given path
 * to accurate evaluation, it uses
 * Simpson's one-third rule
 * @param t time
 * @param zeta path (zeta part)
 * @param theta path (theta part)
 * @return the action of given path
 */

double eval_action(std::vector<double> &t,
fourier_path &zeta, fourier_path &theta);
