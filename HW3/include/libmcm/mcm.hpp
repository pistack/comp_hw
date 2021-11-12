/*!
 * @file mcm.hpp
 * @ingroup libmcm
 * @brief headerfile for the
 * minimization of the action by Monte Carlo Metropolis
 * method described in 
 * [Entropy 2020, 22(9), 916](https://doi.org/10.3390/e22090916)
 * @author pistack (Junho Lee)
 * @date 2021. 11. 12.
 */

#ifndef MCM_H
#define MCM_H

#include <cmath>
#include <cstddef>
#include <fstream>
#include <tuple>
#include <vector>
#include <string>
#include <random>
#include "libpath/action.hpp"

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
/// @param Basis type of basis used to approximate
/// path
/// @param Path type of path
/// @param Lag lagranian of action
/// functor class which has
/// time, path and derivative of path
/// as variable and it returns
/// value of lagranian at given time
/// @see \ref mcm
/// @note This is abstract class use
/// libmcm::mcm_fourier (for libpath::fourier_path) or
/// libmcm::mcm_bezier (for libpath::bezier_path), instead.
/// @ingroup libmcm

template<typename T, typename Basis, typename Path, typename Lag>
class mcm
{
   private:

   T t0, t1; // initial and final time
   std::vector<T> p0, p1; // initial and final value of path

   // setup for basis
   int order;
   T add_setup;

   // action
   libpath::action<T, Path, Lag> mcm_action;

   // initial guess
   std::vector<std::vector<T>> init_guess;
   std::vector<Path> init_path;

   // result
   std::vector<Path> min_path;

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
   rand_walk(std::vector<std::vector<T>> guess,
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
   /// @param order_ order of basis function
   mcm(T t_0, T t_1, 
   std::vector<T> p_0, std::vector<T> p_1,
   T abs_tol, int order_)
   : t0(t_0), t1(t_1), p0(p_0), p1(p_1), \
   order(order_), add_setup(0)
   {
      mcm_action = libpath::action<T, Path, Lag>(abs_tol);
   }

   /// @brief initialize mcm class
   /// @param t_0 initial time
   /// @param t_1 finial time
   /// @param p_0 value of path at initial time
   /// @param p_1 value of path at finial time
   /// @param abs_tol absolute tolerance for action integral
   /// @param order_ order of basis function
   /// @param add_setup_ additional setup used to define basis function
   mcm(T t_0, T t_1, 
   std::vector<T> p_0, std::vector<T> p_1,
   T abs_tol, int order_, T add_setup_)
   : t0(t_0), t1(t_1), p0(p_0), p1(p_1), \
   order(order_), add_setup(add_setup_)
   {
      mcm_action = libpath::action<T, Path, Lag>(abs_tol);
   }

   /// @brief copy constructer of mcm class
   mcm(const mcm<T, Basis, Path, Lag> &copy)
   : t0(copy.t0), t1(copy.t1), p0(copy.p0), p1(copy.p1), \
   order(copy.order), \
   add_setup(copy.add_setup), \
   mcm_action(copy.mcm_action), \
   init_path(copy.init_path), \
   min_path(copy.min_path)
   {}

   /// @brief overloading of assignment operator for 
	/// mcm class
   mcm<T, Basis, Path, Lag> & operator=(const mcm<T, Basis, Path, Lag> &copy);

   /// @brief set initial guess
   /// @param init_c initial coefficients
   void set_init_guess(std::vector<std::vector<T>> init_c);

   /// @brief set initial guess randomly
   void set_init_guess();

   /// @brief get action of initial guess
   /// @param[out] e estimated error of action integral
   /// @return action of initial guess
   /// @see libpath::action 
   T get_init_action(T &e);

   /// @brief get initial path
   /// @return initial_path
   std::vector<Path> get_init_path() const
   {return init_path;}

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
   /// @see libpath::action
   T get_min_action(T &e);

   /// @brief get minimum path
   /// @return minimum path
   std::vector<Path>
   get_min_coeff() const
   {return min_path;}

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
   std::tuple<std::size_t, T>
   optimize(std::size_t num_iter, T step_size, T lambda);

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
   std::tuple<std::size_t, T>
   optimize(std::size_t num_iter, T step_size, T lambda, 
   std::string monitor);
};
}

#include "mcm/mcm_basic.tpp"
#include "mcm/mcm_move.tpp"
#include "mcm/mcm_opt.tpp"

#endif