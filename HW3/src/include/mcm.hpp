/*!
 * @file mcm.hpp
 * @ingroup libmcm
 * @brief headerfile to
 * minimize the action by Monte Carlo Metropolis
 * method described in 
 * [Entropy 2020, 22(9), 916](https://doi.org/10.3390/e22090916)
 * @author pistack (Junho Lee)
 * @date 2021. 10. 26.
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

/// @brief class to
/// Minimize the action via Monte Carlo
/// Metropolis Method.
/// It uses mt19937 random number generator to
/// generates distribution.
/// Moreover it samples path via random walk.
/// To make random walk, it samples real number
/// from normal distribution and move path by
/// the sampled real number.
/// If type T is not one of float, double and
/// long doulbe then class is undefined.
/// @see \ref mcm
/// @ingroup libmcm

template<typename T>
class mcm
{
   private:

   T t0, t1; // initial and final time
   std::vector<T> p0, p1; // initial and final value of path

   // setup for fourier function
   int num_fourier;
   T fourier_period;

   // action
   action<T> mcm_action;

   // initial guess
   std::vector<std::vector<T>> init_guess;
   std::vector<fourier_path<T>> init_path;
   T init_action;

   // result
   std::vector<std::vector<T>> min_guess;
   std::vector<fourier_path<T>> min_path;
   T min_action;

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
   /// @param rel_tol relative tolerance for action integral
   /// @param num_f number of sine and cosine function to use
   /// @param period period of fourier function
   /// @param lag lagranian of action
   /// @see fourier class and action class
   mcm(T t_0, T t_1, 
   std::vector<T> p_0, std::vector<T> p_1,
   T abs_tol, T rel_tol, int num_f, T period,
   T (*lag)(T, std::vector<T>, 
   std::vector<T>))
   : t0(t_0), t1(t_1), p0(p_0), p1(p_1), \
   num_fourier(num_f), fourier_period(period)
   {
      mcm_action = action<T>(abs_tol, rel_tol, lag);
   }

   /// @brief copy constructer of HW3 class
   mcm(const mcm<T> &copy)
   : t0(copy.t0), t1(copy.t1), p0(copy.p0), p1(copy.p1), \
   num_fourier(copy.num_fourier), \
   fourier_period(copy.fourier_period), \
   mcm_action(copy.mcm_action), init_guess(copy.init_guess), \
   init_path(copy.init_path), init_action(copy.init_action), \
   min_guess(copy.min_guess), min_path(copy.min_path), \
   min_action(copy.min_action)
   {}

   mcm<T> & operator=(const mcm<T> &copy);

   /// @brief set initial guess
   /// @param init_c initial coefficients to weight sum of
   /// sine and cosine function
   /// @see fourier class
   void set_init_guess(std::vector<std::vector<T>> init_c);

   /// @brief set initial guess randomly
   void set_init_guess();

   /// @brief get action of initial guess
   /// @return action of initial guess
   T get_init_action();

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
   /// @return action of minimum guess
   T get_min_action();

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
   /// @return tuple of number of accepted move and acceptance ratio.
   std::tuple<int, T>
   optimize(int num_iter, T step_size, T lambda, 
   std::string monitor);
};

#include "mcm/mcm_basic.tpp"
#include "mcm/mcm_move.tpp"
#include "mcm/mcm_opt.tpp"

#endif