/*!
 * @file test_action_simple.cpp
 * @brief test action::eval() routine with simple function
 * @author pistack (Junho Lee)
 * @date 2021. 10. 30.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include "action.hpp"

#if PRECISION_LEVEL == 0
    #define PRECISION float
    #define DIGITS 6
#elif PRECISION_LEVEL == 1
    #define PRECISION double
    #define DIGITS 14
#endif

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
class id_inv_lag{
  public :
  id_inv_lag() {}
  T operator()(T t, vector<T> p, vector<T> dp) const
  {
    return p[0]+1/p[0];
  }
};

int main(void)
{

  cout << "==========================================================" << endl;
  cout << "               Test eval action routine                   " << endl;
  cout << " Test 1. integrate 2.0 + sin(pi*x)+cos(pi*x) from 0 to 2                " << endl;
  cout << " Test 2. integrate 1/(2.0+sin(pi*x)+cos(pi*x)) from 0 to 2    " << endl;
  cout << " Test 3. Test 1 + Test 2 = 2.0 + Test 2 " << endl;

  // initial condition
  vector<PRECISION> c = {1.0, 1.0};
  fourier<PRECISION> tmp(1, 2.0, c);
  vector<fourier_path<PRECISION>> path(1, fourier_path<PRECISION>(0.0, 2.0, 3.0, 3.0, tmp));
  action<PRECISION, id_lag<PRECISION>> id_action1(1e-4, path);
  action<PRECISION, id_lag<PRECISION>> id_action2(1e-6, path);
  action<PRECISION, id_lag<PRECISION>> id_action3(1e-8, path);
  action<PRECISION, inv_lag<PRECISION>> inv_action1(1e-4, path);
  action<PRECISION, inv_lag<PRECISION>> inv_action2(1e-6, path);
  action<PRECISION, inv_lag<PRECISION>> inv_action3(1e-8, path);
  action<PRECISION, id_inv_lag<PRECISION>> id_inv_action1(1e-4, path);
  action<PRECISION, id_inv_lag<PRECISION>> id_inv_action2(1e-6, path);
  action<PRECISION, id_inv_lag<PRECISION>> id_inv_action3(1e-8, path);

  cout.unsetf(ios::floatfield); // initialize floatfield
  cout.precision(DIGITS); // print significant digits
  cout << " Test 1. atol: " << 1e-4 <<  endl;
  cout << "Integration value: " << id_action1.eval() << endl;
  cout << " Test 1. atol: " << 1e-6 << endl;
  cout << "Integration value: " << id_action2.eval() << endl;
  cout << " Test 1. atol: " << 1e-8 << endl;
  cout << "Integration value: " << id_action3.eval() << endl;
  cout << " Test 2. atol: " << 1e-4 <<  endl;
  cout << "Integration value: " << inv_action1.eval() << endl;
  cout << " Test 2. atol: " << 1e-6 << endl;
  cout << "Integration value: " << inv_action2.eval() << endl;
  cout << " Test 2. atol: " << 1e-8 << endl;
  cout << "Integration value: " << inv_action3.eval() << endl;
  cout << " Test 3. atol: " << 1e-4 <<  endl;
  cout << "Integration value: " << id_inv_action1.eval() << endl;
  cout << " Test 3. atol: " << 1e-6 << endl;
  cout << "Integration value: " << id_inv_action2.eval() << endl;
  cout << " Test 3. atol: " << 1e-8 << endl;
  cout << "Integration value: " << id_inv_action3.eval() << endl;
  cout << "Test finished!" << endl;
  cout << "==========================================================" << endl;

  return 0;
}
      
      

  
						  
