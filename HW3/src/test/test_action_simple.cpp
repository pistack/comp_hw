/*!
 * @file test_action_simple.cpp
 * @brief test action::eval() routine with simple function
 * @author pistack (Junho Lee)
 * @date 2021. 11. 2.
 */

#include <chrono>
#include <cmath>
#include <iostream>
#include "action.hpp"

#if PRECISION_LEVEL == 0
    #define PRECISION float
    #define DIGITS 6
#elif PRECISION_LEVEL == 1
    #define PRECISION double
    #define DIGITS 14
#endif

using namespace libfourier;
using namespace std;

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
class id_lag{
  public :
  id_lag() {}
  T operator()(T t, vector<T> p, vector<T> dp) const
  {
    return p[0];
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
class square_lag{
  public :
  square_lag() {}
  T operator()(T t, vector<T> p, vector<T> dp) const
  {
    return std::pow(p[0], 2.0);
  }
};

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
class weird_lag{
  public :
  weird_lag() {}
  T operator()(T t, vector<T> p, vector<T> dp) const
  {
    return std::sin(2*std::exp(2*(std::sin(2*(std::exp(2*p[0]))))));
  }
};

int main(void)
{

  cout << "==========================================================" << endl;
  cout << "               Test action::eval() routine                " << endl;
  cout << " Test 1. integrate 2.0 + sin(pi*x) from 0 to 2                " << endl;
  cout << " Test 2. integrate 1/(2.0+sin(pi*x)) from 0 to 2    " << endl;
  cout << " Test 3. integrate sin(sin(pi*x)) from 0 to 2 " << endl;
  cout << " Test 4. integrate sin(pi*x)**2 from 0 to 2 " << endl;
  cout << " Test 5. integrate 0 from 0 to 2 " << endl;
  cout << " Test 6. integrate sin(2*exp(2*sin(2(exp(2*sin(pi*x)))))) from 0 to 2 " << endl;

  // order of Gauss-Kronrod quadrature method
  vector<int> order = {15, 21, 31, 41, 51, 61};
  // execution time
  std::chrono::steady_clock::time_point start;
  std::chrono::steady_clock::time_point end;

  // initial condition
  vector<PRECISION> c = {1.0, 0.0};
  #if PRECISION_LEVEL == 0
  vector<PRECISION> tol = {1.0, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10};
  #elif PRECISION_LEVEL == 1
  vector<PRECISION> tol = {1.0, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14};
  #endif
  fourier<PRECISION> tmp(1, 2.0, c);
  vector<fourier_path<PRECISION>> path(1, fourier_path<PRECISION>(0.0, 2.0, 2.0, 2.0, tmp));
  vector<fourier_path<PRECISION>> path2(1, fourier_path<PRECISION>(0.0, 2.0, 0.0, 0.0, tmp));
  action<PRECISION, id_lag<PRECISION>> tst1(path);
  action<PRECISION, inv_lag<PRECISION>> tst2(path);
  action<PRECISION, sine_lag<PRECISION>> tst3(path2);
  action<PRECISION, square_lag<PRECISION>> tst4(path2);
  action<PRECISION, zero_lag<PRECISION>> tst5(path2);
  action<PRECISION, weird_lag<PRECISION>> tst6(path2);
  cout.unsetf(ios::floatfield); // initialize floatfield
  cout.precision(DIGITS); // print significant digits
  for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst1.update(*it);
    cout << " Test 1. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst1.eval();
    end = std::chrono::steady_clock::now();
    cout << "Integration value: " << tst1.eval() << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst2.update(*it);
    cout << " Test 2. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst2.eval();
    end = std::chrono::steady_clock::now();
    cout << "Integration value: " << tst2.eval() << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst3.update(*it);
    cout << " Test 3. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<10000; ++j)
    tst3.eval();
    end = std::chrono::steady_clock::now();
    cout << "Integration value: " << tst3.eval() << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst4.update(*it);
    cout << " Test 4. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<1; ++j)
    tst4.eval();
    end = std::chrono::steady_clock::now();
    cout << "Integration value: " << tst4.eval() << '\t' << tst4.eval(1, 10) << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst5.update(*it);
    cout << " Test 5. atol: " << *it <<  endl;
    start = std::chrono::steady_clock::now();
    for(int j=0; j<1; ++j)
    tst5.eval();
    end = std::chrono::steady_clock::now();
    cout << "Integration value: " << tst5.eval() << '\t' << tst5.eval(1, 10) << endl;
    cout << "Execution time: " << \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/10000.0 << \
    " microsecond" << endl;
  }
  for(std::vector<PRECISION>::iterator it=tol.begin(); it != tol.end(); ++it)
  {
    tst6.update(*it);
    cout << " Test 6. atol: " << *it <<  endl;
    cout << "Integration Method: Gauss-Kron quadrature" << endl;
    for(int i=0; i<4; ++i)
    {
      start = std::chrono::steady_clock::now();
      for(int j=0; j<1000; ++j)
      tst6.eval(order[i]);
      end = std::chrono::steady_clock::now();
      cout << "Order: " << order[i] << endl;
      cout << "Integration value: " << tst6.eval(order[i]) << endl;
      cout << "Execution time: " << \
      std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000.0 << \
      " microsecond" << endl;
    }
  }
  cout << "Test finished!" << endl;
  cout << "==========================================================" << endl;

  return 0;
}
      
      

  
						  
