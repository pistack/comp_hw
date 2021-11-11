/*!
 * @file test_action_simple.cpp
 * @brief test action::eval() routine with simple function
 * @author pistack (Junho Lee)
 * @date 2021. 11. 11.
 */

#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include "libpath/bezier_path.hpp"
#include "libpath/fourier_path.hpp"
#include "libpath/action.hpp"

#if PRECISION_LEVEL == 0
    #define PRECISION float
    #define DIGITS 6
#elif PRECISION_LEVEL == 1
    #define PRECISION double
    #define DIGITS 14
#endif

using namespace libpath;
using namespace std;

#ifndef DOXYGEN_SKIP

template<typename T>
class zero_lag{
  public :
  zero_lag() {}
  T operator()(T t, vector<T> p, vector<T> dp) const
  {
    return 0;
  }
};

template<typename T>
class id_lag{
  public :
  id_lag() {}
  T operator()(T t, vector<T> p, vector<T> dp) const
  {
    return p[0];
  }
};

template<typename T>
class inv_lag{

  public :

  inv_lag() {}

  T operator()(T t, 
  vector<T> p, vector<T> dp) const
  {
    return 1/p[0];
  }
};

template<typename T>
class square_lag{
  public :
  square_lag() {}
  T operator()(T t, vector<T> p, vector<T> dp) const
  {
    return p[0]*p[0];
  }
};

template<typename T>
class sine_lag{
  public :
  sine_lag() {}
  T operator()(T t, vector<T> p, vector<T> dp) const
  {
    return std::sin(p[0]);
  }
};

template<typename T>
class exp_lag{
  public :
  exp_lag() {}
  T operator()(T t, vector<T> p, vector<T> dp) const
  {
    return std::exp(-p[0]);
  }
};

template<typename T>
class inv_sqrt_lag{

  public :

  inv_sqrt_lag() {}

  T operator()(T t, 
  vector<T> p, vector<T> dp) const
  {
    return 1/std::sqrt(p[0]);
  }
};

template<typename T>
class weird_lag{
  public :
  weird_lag() {}
  T operator()(T t, vector<T> p, vector<T> dp) const
  {
    return std::sin(2*std::exp(2*(std::sin(2*(std::exp(2*p[0]))))));
  }
};
#endif /* DOXYGEN_SKIP */

int main(void)
{

  cout << "==============================================================" << endl;
  cout << "       Test action::eval() routine simple lagrangian          " << endl;
  cout << " Test 0. integrate 0 from 0 to 2                              " << endl;
  cout << " Test 1. integrate sin(pi*x) from 0 to 1                      " << endl;
  cout << " Test 2. integrate sin(pi*x) from 0 to 2                      " << endl;
  cout << " Test 3. integrate 1/(2.0+sin(pi*x)) from 0 to 1              " << endl;
  cout << " Test 4. integrate 1/(2.0+sin(pi*x)) from 0 to 2              " << endl;
  cout << " Test 5. integrate sin(pi*x)**2 from 0 to 1                   " << endl;
  cout << " Test 6. integrate sin(sin(pi*x)) from 0 to 2                 " << endl;
  cout << " Test 7. integrate exp(sin(pi*x)) from 0 to 2                 " << endl;
  cout << " Test 8. integrate x from 0 to 1                              " << endl;
  cout << " Test 9. integrate 10th order bezier curve from 0 to 1        " << endl;
  cout << " Control points: (1,2,3,4,5,6,5,4,3,2,1)                      " << endl;
  cout << " Test 10. integrate 1/(1+x) from 0 to 1                       " << endl;
  cout << " Test 11. integrate 1/(1+x^2) from 0 to 1                     " << endl;
  cout << " Test 12. integrate exp(-x**2) from 0 to 1                    " << endl;
  cout << " Test Integrand which has singularties at end points          " << endl;
  cout << " Test 13. integrate 1/sqrt(sin(pi*x)) from 0 to 1             " << endl;
  cout << " Test 14. integrate 1/sqrt(x) from 0 to 1                     " << endl;
  cout << " Test 15. integrate 1/sqrt((x-0.9)*(9/8-x)) from 0.9 to 9/8   " << endl;
  cout << " Test wildly oscillate integrand                              " << endl;
  cout << " Test 16. integrate sin(2*exp(2*sin(2*exp(2*sin(pi*x))))) from 0 to 2 " << endl;

  // exact result
  // values are obtained by wolfram-alpha
  vector<PRECISION> exact {
    0, 2/pi<PRECISION>, 0, 0.63661977236758134308, 1.15470053837925152902, 0.5, 0, 
    2.5321317555040166712,
    0.5, 36.0/11, std::log(2.0), pi<PRECISION>/4, std::sqrt(pi<PRECISION>)*std::erf(1.0)/2,
    1.6692536833481463726, 2, pi<PRECISION>,
    0.0561899582642922203122
  };

  // order of Gauss-Kronrod quadrature method
  vector<int> order = {15, 21, 31, 41, 51, 61};
  // execution time
  std::chrono::steady_clock::time_point start;
  std::chrono::steady_clock::time_point end;
  // estimated error
  PRECISION e;

  // initial condition
  // fourier func
  PRECISION T = 2; // period
  vector<PRECISION> c = {1, 0}; // use one sine function
  // bezier curve
  vector<PRECISION> B_c1 = {0, 1}; // x
  vector<PRECISION> B_c2 = {1, 2}; // 1+x
  vector<PRECISION> B_c3 = {1, 2, 2}; // 1+x^2
  vector<PRECISION> B_c4 = {0, 0, 1}; // x^2
  vector<PRECISION> B_c5 = {0, 1, 0}; // x(1-x)
  vector<PRECISION> B_c6 = {1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1}; // 10th order bezier

  #if PRECISION_LEVEL == 0
  vector<PRECISION> tol = {1.0, 1e-2, 1e-4};
  #elif PRECISION_LEVEL == 1
  vector<PRECISION> tol = {1.0, 1e-2, 1e-4, 1e-6, 1e-8};
  #endif
  // bezier
  bezier<PRECISION> B_fun1(1, B_c1); // x
  bezier<PRECISION> B_fun2(1, B_c2); // 1+x
  bezier<PRECISION> B_fun3(2, B_c3); // 1+x^2
  bezier<PRECISION> B_fun4(1, B_c4); // x^2
  bezier<PRECISION> B_fun5(1, B_c5); // x(1-x)
  bezier<PRECISION> B_fun6(1, B_c6); // 10 th order bezier
  // x
  vector<bezier_path<PRECISION>> B_path1(1, bezier_path<PRECISION>(0, 1, 0, 1, B_fun1));
  // 1+x
  vector<bezier_path<PRECISION>> B_path2(1, bezier_path<PRECISION>(0, 1, 1, 2, B_fun2));
  // 1+x^2
  vector<bezier_path<PRECISION>> B_path3(1, bezier_path<PRECISION>(0, 1, 1, 2, B_fun3));
  // -x^2
  vector<bezier_path<PRECISION>> B_path4(1, bezier_path<PRECISION>(0, 1, 0, -1, B_fun4));
  // (x-0.9)(9/8-x)
  vector<bezier_path<PRECISION>> B_path5(1, bezier_path<PRECISION>(0.9, 9./8, 0, 0, B_fun5));
  // 10 th order bezier
  vector<bezier_path<PRECISION>> B_path6(1, bezier_path<PRECISION>(0, 1, 1, 1, B_fun6));
  // fourier
  fourier<PRECISION> tmp(1, T, c);
  // sin(pi*x) [0, 1]
  vector<fourier_path<PRECISION>> f_path1(1, fourier_path<PRECISION>(0, 1, 0, 0, tmp));
  // sin(pi*x) [0, 2]
  vector<fourier_path<PRECISION>> f_path2(1, fourier_path<PRECISION>(0, 2, 0.0, 0.0, tmp));
  // 2+sin(pi*x) [0, 1]
  vector<fourier_path<PRECISION>> f_path3(1, fourier_path<PRECISION>(0, 1, 2, 2, tmp));
  // 2+sin(pi*x) [0, 2]
  vector<fourier_path<PRECISION>> f_path4(1, fourier_path<PRECISION>(0, 2, 2, 2, tmp));

  // test cases
  // simple fourier
  action<PRECISION, fourier_path<PRECISION>, zero_lag<PRECISION>> tst0(f_path2);
  action<PRECISION, fourier_path<PRECISION>, id_lag<PRECISION>> tst1(f_path1);
  action<PRECISION, fourier_path<PRECISION>, id_lag<PRECISION>> tst2(f_path2);
  action<PRECISION, fourier_path<PRECISION>, inv_lag<PRECISION>> tst3(f_path3);
  action<PRECISION, fourier_path<PRECISION>, inv_lag<PRECISION>> tst4(f_path4);
  action<PRECISION, fourier_path<PRECISION>, square_lag<PRECISION>> tst5(f_path1);
  action<PRECISION, fourier_path<PRECISION>, sine_lag<PRECISION>> tst6(f_path2);
  action<PRECISION, fourier_path<PRECISION>, exp_lag<PRECISION>> tst7(f_path2);
  // simple bezier
  action<PRECISION, bezier_path<PRECISION>, id_lag<PRECISION>> tst8(B_path1);
  action<PRECISION, bezier_path<PRECISION>, id_lag<PRECISION>> tst9(B_path6);
  action<PRECISION, bezier_path<PRECISION>, inv_lag<PRECISION>> tst10(B_path2);
  action<PRECISION, bezier_path<PRECISION>, inv_lag<PRECISION>> tst11(B_path3);
  action<PRECISION, bezier_path<PRECISION>, exp_lag<PRECISION>> tst12(B_path4);

  // singularity
  action<PRECISION, fourier_path<PRECISION>, inv_sqrt_lag<PRECISION>> tst13(f_path1);
  action<PRECISION, bezier_path<PRECISION>, inv_sqrt_lag<PRECISION>> tst14(B_path1);
  action<PRECISION, bezier_path<PRECISION>, inv_sqrt_lag<PRECISION>> tst15(B_path5);

  // performance test
  action<PRECISION, fourier_path<PRECISION>, weird_lag<PRECISION>> tst16(f_path2);

  // tst setup
  PRECISION result;
  PRECISION exact_error;
  bool sucess = true;
  int tst_num = 0;
  
  cout.unsetf(ios::floatfield); // initialize floatfield
  cout.precision(DIGITS); // print significant digits
  // simple test
  // test0
  for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst0.update(*it);
    cout << " Test 0. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst0.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: default (G15, K31) gauss-kronrod quadrature" << endl;
    result = tst0.eval(e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    for(int j=0; j<10000; ++j)
    tst0.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: Tanh-Sinh Quadrature" << endl;
    result = tst0.eval(1, 7, e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  // test 1
    for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst1.update(*it);
    cout << " Test 0. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: default (G15, K31) gauss-kronrod quadrature" << endl;
    result = tst1.eval(e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    for(int j=0; j<10000; ++j)
    tst1.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: Tanh-Sinh Quadrature" << endl;
    result = tst1.eval(1, 7, e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  ++tst_num;
  // test 2
    for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst1.update(*it);
    cout << " Test 0. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: default (G15, K31) gauss-kronrod quadrature" << endl;
    result = tst1.eval(e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    for(int j=0; j<10000; ++j)
    tst1.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: Tanh-Sinh Quadrature" << endl;
    result = tst1.eval(1, 7, e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  ++tst_num;
  // test 3
    for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst1.update(*it);
    cout << " Test 0. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: default (G15, K31) gauss-kronrod quadrature" << endl;
    result = tst1.eval(e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    for(int j=0; j<10000; ++j)
    tst1.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: Tanh-Sinh Quadrature" << endl;
    result = tst1.eval(1, 7, e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  ++tst_num;
  // test 4
    for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst1.update(*it);
    cout << " Test 0. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: default (G15, K31) gauss-kronrod quadrature" << endl;
    result = tst1.eval(e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    for(int j=0; j<10000; ++j)
    tst1.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: Tanh-Sinh Quadrature" << endl;
    result = tst1.eval(1, 7, e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  ++tst_num;
  // test 5
    for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst1.update(*it);
    cout << " Test 0. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: default (G15, K31) gauss-kronrod quadrature" << endl;
    result = tst1.eval(e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    for(int j=0; j<10000; ++j)
    tst1.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: Tanh-Sinh Quadrature" << endl;
    result = tst1.eval(1, 7, e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  ++tst_num;
  // test 6
    for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst1.update(*it);
    cout << " Test 0. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: default (G15, K31) gauss-kronrod quadrature" << endl;
    result = tst1.eval(e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    for(int j=0; j<10000; ++j)
    tst1.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: Tanh-Sinh Quadrature" << endl;
    result = tst1.eval(1, 7, e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  ++tst_num;
  // test 7
    for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst1.update(*it);
    cout << " Test 0. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: default (G15, K31) gauss-kronrod quadrature" << endl;
    result = tst1.eval(e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    for(int j=0; j<10000; ++j)
    tst1.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: Tanh-Sinh Quadrature" << endl;
    result = tst1.eval(1, 7, e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  ++tst_num;
  // test 8
    for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst1.update(*it);
    cout << " Test 0. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: default (G15, K31) gauss-kronrod quadrature" << endl;
    result = tst1.eval(e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    for(int j=0; j<10000; ++j)
    tst1.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: Tanh-Sinh Quadrature" << endl;
    result = tst1.eval(1, 7, e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  ++tst_num;
  // test 9 
    for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst1.update(*it);
    cout << " Test 0. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: default (G15, K31) gauss-kronrod quadrature" << endl;
    result = tst1.eval(e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    for(int j=0; j<10000; ++j)
    tst1.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: Tanh-Sinh Quadrature" << endl;
    result = tst1.eval(1, 7, e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  ++tst_num;
  // test 10
    for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst1.update(*it);
    cout << " Test 0. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: default (G15, K31) gauss-kronrod quadrature" << endl;
    result = tst1.eval(e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    for(int j=0; j<10000; ++j)
    tst1.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: Tanh-Sinh Quadrature" << endl;
    result = tst1.eval(1, 7, e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  ++tst_num;
  // test 11
    for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst1.update(*it);
    cout << " Test 0. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: default (G15, K31) gauss-kronrod quadrature" << endl;
    result = tst1.eval(e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    for(int j=0; j<10000; ++j)
    tst1.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: Tanh-Sinh Quadrature" << endl;
    result = tst1.eval(1, 7, e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  ++tst_num;
  // test 12
    for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst1.update(*it);
    cout << " Test 0. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval(e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: default (G15, K31) gauss-kronrod quadrature" << endl;
    result = tst1.eval(e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
    for(int j=0; j<10000; ++j)
    tst1.eval(1, 7, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: Tanh-Sinh Quadrature" << endl;
    result = tst1.eval(1, 7, e);
    exact_error = std::abs(result-exact[tst_num]);
    if(exact_error>(*it))
    sucess = false;
    cout << "Integration value: " << result << endl;
    cout << "Exact value from Wolfram-Alpha: " << exact[tst_num] << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Exact error: " << exact_error << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  ++tst_num;
  // singularity in end points
  for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst7.update(*it);
    cout << " Test 7. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<1000; ++j)
    tst7.eval(1, 10, e);
    end = std::chrono::steady_clock::now();
    cout << "Integration Method: Tanh-Sinh quadrature" << endl;
    cout << "Integration value: " << tst7.eval(1, 10, e) << endl;
    cout << "Estimated error: " << e << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000.0 << \
    " microsecond" << endl;
  }

  // perforcemance test
  for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst6.update(*it);
    cout << " Test 6. atol: " << *it <<  endl;
    cout << "Integration Method: Gauss-Kron quadrature" << endl;
    for(int i=0; i<6; ++i)
    {
      start = std::chrono::steady_clock::now();
      for(int j=0; j<1000; ++j)
      tst6.eval(0, order[i], e);
      end = std::chrono::steady_clock::now();
      cout << "Order: " << order[i] << endl;
      cout << "Integration value: " << tst6.eval(0, order[i], e) << endl;
      cout << "Estimated error: " << e << endl;
      cout << "Execution time: " << \
      std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000.0 << \
      " microsecond" << endl;
    }
  }
  cout << "Test finished!" << endl;
  cout << "==========================================================" << endl;

  return 0;
}
      
      

  
						  
