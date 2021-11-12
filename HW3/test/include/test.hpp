/*!
 * @file test.hpp
 * @brief header file for testing libpath
 * @author pistack (Junho Lee)
 * @date 2021. 11. 12.
 */

#include <vector>
#include <cmath>

#ifndef TEST_H
#define TEST_H


#if PRECISION_LEVEL == 0
    #define PRECISION float
    #define DIGITS 6
#elif PRECISION_LEVEL == 1
    #define PRECISION double
    #define DIGITS 14
#endif

namespace test {

#ifndef DOXYGEN_SKIP
template<typename T>
class kepler_lag{

  public :

  kepler_lag() {}

  T operator()(T t, 
  std::vector<T> p, std::vector<T> dp) const
  {
    return 0.5*(std::pow(dp[0], 2.0)+std::pow(p[0]*dp[1], 2.0))+
    1/std::abs(p[0]);
  }

};

template<typename T>
class zero_lag{
  public :
  zero_lag() {}
  T operator()(T t, std::vector<T> p, std::vector<T> dp) const
  {
    return 0;
  }
};

template<typename T>
class id_lag{
  public :
  id_lag() {}
  T operator()(T t, std::vector<T> p, std::vector<T> dp) const
  {
    return p[0];
  }
};

template<typename T>
class inv_lag{

  public :

  inv_lag() {}

  T operator()(T t, 
  std::vector<T> p, std::vector<T> dp) const
  {
    return 1/p[0];
  }
};

template<typename T>
class square_lag{
  public :
  square_lag() {}
  T operator()(T t, std::vector<T> p, std::vector<T> dp) const
  {
    return p[0]*p[0];
  }
};

template<typename T>
class sine_lag{
  public :
  sine_lag() {}
  T operator()(T t, std::vector<T> p, std::vector<T> dp) const
  {
    return std::sin(p[0]);
  }
};

template<typename T>
class exp_lag{
  public :
  exp_lag() {}
  T operator()(T t, std::vector<T> p, std::vector<T> dp) const
  {
    return std::exp(p[0]);
  }
};

template<typename T>
class inv_sqrt_lag{

  public :

  inv_sqrt_lag() {}

  T operator()(T t, 
  std::vector<T> p, std::vector<T> dp) const
  {
    return 1/std::sqrt(p[0]);
  }
};

template<typename T>
class weird_lag{
  public :
  weird_lag() {}
  T operator()(T t, std::vector<T> p, std::vector<T> dp) const
  {
    return std::sin(2*std::exp(2*(std::sin(2*(std::exp(2*p[0]))))));
  }
};
#endif

/// @brief test libpath::action::eval() routine with kepler lagrangian
/// @return test result
int test_action_kepler();

/// @brief test libpath::action::eval() routine with various simple lagrangian
/// @return test result
int test_action_simple();

/// @brief test libpath::action::is_vaild() routine 
/// @return test result
int test_action_vaildity();

/// @brief test libpath::bezier_path class 
/// @return test result
int test_bezier_path();

/// @brief test libpath::bezier class 
/// @return test result
int test_bezier();

/// @brief test libpath::fourier_path class
/// @return test result
int test_fourier_path();

/// @brief test libpath::fourier class
/// @return test result
int test_fourier();
}
 
#endif