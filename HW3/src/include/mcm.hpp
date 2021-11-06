/*!
 * @file mcm.hpp
 * @ingroup libmcm
 * @brief headerfile for the
 * minimization of the action by Monte Carlo Metropolis
 * method described in 
 * [Entropy 2020, 22(9), 916](https://doi.org/10.3390/e22090916)
 * @author pistack (Junho Lee)
 * @date 2021. 11. 7.
 */

#ifndef MCM_H
#define MCM_H

#include <cmath>
#include <fstream>
#include <tuple>
#include <vector>
#include <string>
#include <random>
#include "fourier_path.hpp"
#include "action.hpp"

namespace libmcm {
/// @brief class to
/// Minimize the action via Monte Carlo
/// Metropolis Method.
/// It uses mt19937 random number generator to
/// generates distribution.
/// Moreover it samples path via random walk.
/// To make random walk, it samples real number
/// from normal distribution and move path by
/// the sampled real number.
/// @param T precision should be one of
/// float, double and long double
/// @param Lag lagranian of action
/// functor class which has
/// time, path and derivative of path
/// as variable and it returns
/// value of lagranian at given time
/// @see \ref mcm
/// @ingroup libmcm

template<typename T, typename Lag>
class mcm
{
   private:

   T t0, t1; // initial and final time
   std::vector<T> p0, p1; // initial and final value of path

   // setup for fourier function
   int num_fourier;
   T fourier_period;

   // action
   libfourier::action<T, Lag> mcm_action;

   // initial guess
   std::vector<std::vector<T>> init_guess;
   std::vector<libfourier::fourier_path<T>> init_path;

   // result
   std::vector<std::vector<T>> min_guess;
   std::vector<libfourier::fourier_path<T>> min_path;

   std::random_device rd;
   std::mt19937 gen = std::mt19937(rd()); // set random number generator
   std::normal_distribution<T> normal_dist = \
   std::normal_distribution<T>(0.0, 1.0); // distribution for random walk
   std::uniform_real_distribution<T> uniform_dist = \
   std::uniform_real_distribution<T>(0.0, 1.0); // set distribution for random number

   /// @brief sample guess via random walk
   /// @param guess guess to move 
   /// @param step_size step_size of random walk
   /// @return moved guess via random walk
   std::vector<std::vector<T>>
   move(std::vector<std::vector<T>> guess,
   T step_size);

   public:

   /// @brief initialize mcm class
   mcm(){};

   /// @brief initialize mcm class
   /// @param t_0 initial time
   /// @param t_1 finial time
   /// @param p_0 value of path at initial time
   /// @param p_1 value of path at finial time
   /// @param abs_tol absolute tolerance for action integral
   /// @param num_f number of sine and cosine function to use
   /// @param period period of fourier function
   /// @see fourier class and action class
   mcm(T t_0, T t_1, 
   std::vector<T> p_0, std::vector<T> p_1,
   T abs_tol, int num_f, T period)
   : t0(t_0), t1(t_1), p0(p_0), p1(p_1), \
   num_fourier(num_f), fourier_period(period)
   {
      mcm_action = libfourier::action<T, Lag>(abs_tol);
   }

   /// @brief copy constructer of HW3 class
   mcm(const mcm<T, Lag> &copy)
   : t0(copy.t0), t1(copy.t1), p0(copy.p0), p1(copy.p1), \
   num_fourier(copy.num_fourier), \
   fourier_period(copy.fourier_period), \
   mcm_action(copy.mcm_action), \
   init_guess(copy.init_guess), init_path(copy.init_path), \
   min_guess(copy.min_guess), min_path(copy.min_path)
   {}

   /// @brief overloading of assignment operator for 
	/// mcm class
   mcm<T, Lag> & operator=(const mcm<T, Lag> &copy);

   /// @brief set initial guess
   /// @param init_c initial coefficients to weight sum of
   /// sine and cosine function
   /// @see libfourier::fourier
   void set_init_guess(std::vector<std::vector<T>> init_c);

   /// @brief set initial guess randomly
   void set_init_guess();

   /// @brief get action of initial guess
   /// @param[out] e estimated error of action integral
   /// @return action of initial guess
   /// @see libfourier::action 
   T get_init_action(T &e);

   /// @brief get coefficients of initial guess
   /// @return tuple of adder, scaler and fourier coefficient
   std::tuple<std::vector<T>, std::vector<T>, 
   std::vector<std::vector<T>>>
   get_init_coeff();

   /// @brief evaluate initial path at given t
   /// @param t time to evaluate initial path
   /// @return initial path evaluated at t
   std::vector<T> init_eval(T t);

   /// @brief evaluate initial path at given t
   /// @param t time to evaluate initial path
   /// @return initial path evaluated at t
   std::vector<std::vector<T>> 
   init_eval(std::vector<T> t);

   /// @brief get action of minimum guess
   /// @param[out] e estimated error of action integral
   /// @return action of minimum guess
   /// @see libfourier::action
   T get_min_action(T &e);

   /// @brief get coefficients of minimum guess
   /// @return tuple of adder, scaler and fourier coefficient  
   std::tuple<std::vector<T>, std::vector<T>, 
   std::vector<std::vector<T>>>
   get_min_coeff();

   /// @brief evaluate minimum path at given t
   /// @param t time to evaluate minimum path
   /// @return minimum path evaluated at t
   std::vector<T> min_eval(T t);

   /// @brief evaluate minium path at given t
   /// @param t time to evaluate minimum path
   /// @return minimum path evaluated at t
   std::vector<std::vector<T>> min_eval(std::vector<T> t);

   /// @brief minimize the action via
   /// Monte Carlo Metropolis method
   /// @param num_iter number of iteration
   /// @param step_size step size of random walk
   /// @param lambda parameter which controls acceptance of move
   /// @return tuple of number of accepted move and acceptance ratio.
   std::tuple<int, T>
   optimize(int num_iter, T step_size, T lambda);

   /// @brief minimize the action via
   /// Monte Carlo Metropolis method
   /// @param num_iter number of iteration
   /// @param step_size step size of random walk
   /// @param lambda parameter which controls acceptance of move
   /// @param monitor filename to monitor optimization process
   /// saved file has c++ binary format.
   /// It stores twice of number of accepted move T type data
   /// as  \f$ a, e, a, e, a, e, \dotsc \f$ ,
   /// where a is the action of \f$ i \f$ th accepted move and
   /// e is the estimated error of \f$ i \f$ th action integration.
   /// @return tuple of number of accepted move and acceptance ratio.
   std::tuple<int, T>
   optimize(int num_iter, T step_size, T lambda, 
   std::string monitor);
};
}

#include "mcm/mcm_basic.tpp"
#include "mcm/mcm_move.tpp"
#include "mcm/mcm_opt.tpp"

#endif