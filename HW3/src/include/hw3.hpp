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
#include "fourier_path.hpp"
#include "action.hpp"

/// @brief class for HW3
/// solve Kepler problem via Markov Chain
/// Monte Carlo Method
/// It uses mt19937 random number generator and
/// uniform distribution form -1 to 1 to
/// sample path.
/// @ingroup hw3

class HW3
{
   private:

   double t0, t1; // initial and final time
   std::vector<double> p0, p1; // initial and final value of path

   // setup for fourier function
   int num_fourier;
   double hw3_period;

   // action
   action hw3_action;

   // initial guess
   std::vector<std::vector<double>> init_guess;
   std::vector<fourier_path> init_path;
   double init_action;

   // result
   std::vector<std::vector<double>> min_guess;
   std::vector<fourier_path> min_path;
   double min_action;

   std::random_device rd;
   std::mt19937 gen = std::mt19937(rd()); // set random number generator
   std::normal_distribution<double> normal_dist = \
   std::normal_distribution<double>(0.0, 1.0); // distribution for random walk
   std::uniform_real_distribution<double> uniform_dist = \
   std::uniform_real_distribution<double>(0.0, 1.0); // set distribution for random number

   /// @brief find the distance of two guess
   /// @param x 
   /// @param y
   /// @return the distance of two guess
   double dist(std::vector<std::vector<double>> x,
   std::vector<std::vector<double>> y);

   /// @brief randomly
   /// move guess by at most max_step
   /// @param guess guess to move 
   /// @param max_step maximum step size to move guess
   /// @return moved guess by at most max_step
   std::vector<std::vector<double>>
   move(std::vector<std::vector<double>> guess,
   double max_step);

   public:

   /// @brief initialize HW3 class
   HW3(){};

   /// @brief initialize HW3 class
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
   HW3(double t_0, double t_1, 
   std::vector<double> p_0, std::vector<double> p_1,
   double abs_tol, double rel_tol, int num_f, double period,
   double (*lag)(double, std::vector<double>, 
   std::vector<double>))
   : t0(t_0), t1(t_1), p0(p_0), p1(p_1), \
   num_fourier(num_f), hw3_period(period)
   {
      hw3_action = action(abs_tol, rel_tol, lag);
   }

   /// @brief copy constructer of HW3 class
   HW3(const HW3 &copy)
   : t0(copy.t0), t1(copy.t1), p0(copy.p0), p1(copy.p1), \
   num_fourier(copy.num_fourier), hw3_period(copy.hw3_period), \
   hw3_action(copy.hw3_action), init_guess(copy.init_guess), \
   init_path(copy.init_path), init_action(copy.init_action), \
   min_guess(copy.min_guess), min_path(copy.min_path), \
   min_action(copy.min_action)
   {}

   HW3 & operator=(const HW3 &copy);

   /// @brief set initial guess
   /// @param init_c initial coefficients to weight sum of
   /// sine and cosine function
   /// @see fourier class
   void set_init_guess(std::vector<std::vector<double>> init_c);

   /// @brief set initial guess randomly
   void set_init_guess();

   /// @brief get action of initial guess
   /// @return action of initial guess
   double get_init_action();

   /// @brief get coefficients of initial guess
   /// @return tuple of adder, scaler and fourier coefficient
   std::tuple<std::vector<double>, std::vector<double>, 
   std::vector<std::vector<double>>>
   get_init_coeff();

   /// @brief evaluate initial path at given t
   /// @param t time to evaluate initial path
   /// @return initial path evaluated at t
   std::vector<double> init_eval(double t);

   /// @brief evaluate initial path at given t
   /// @param t time to evaluate initial path
   /// @return initial path evaluated at t
   std::vector<std::vector<double>> 
   init_eval(std::vector<double> t);

   /// @brief get action of minimum guess
   /// @return action of minimum guess
   double get_min_action();

   /// @brief get coefficients of minimum guess
   /// @return tuple of adder, scaler and fourier coefficient  
   std::tuple<std::vector<double>, std::vector<double>, 
   std::vector<std::vector<double>>>
   get_min_coeff();

   /// @brief evaluate minimum path at given t
   /// @param t time to evaluate minimum path
   /// @return minimum path evaluated at t
   std::vector<double> min_eval(double t);

   /// @brief evaluate minium path at given t
   /// @param t time to evaluate minimum path
   /// @return minimum path evaluated at t
   std::vector<std::vector<double>> min_eval(std::vector<double> t);

   /// @brief optimize path using Markov Chain
   /// Monte Carlo Method
   /// @param max_iter maximum number of iteration
   /// @param max_step maximum step size
   /// @param lambda parameter to accept move
   /// @return tuple of number of accepted move and acceptance ratio.
   std::tuple<int, double>
   optimize(int max_iter, double max_step, double lambda);
};

#endif
