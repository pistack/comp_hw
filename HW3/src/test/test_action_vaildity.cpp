/*!
 * @file test_action_vaildity.cpp
 * @brief test action::is_vaild() routine 
 * @author pistack (Junho Lee)
 * @date 2021. 11. 10.
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include "fourier_path.hpp"
#include "action.hpp"

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
/// dummy lagrangian functor class
/// used to define action
template<typename T>
class zero_lag{
  public :
  zero_lag() {}
  T operator()(T t, vector<T> p, vector<T> dp) const
  {
    return 0;
  }
};
#endif /* DOXYGEN_SKIP */

int main()
{
  bool sucess = true;
  bool result;
  vector<PRECISION> c1 = {1.0, 0.0};
  vector<PRECISION> c2 = {0.0, 1.0};
  vector<fourier_path<PRECISION>> path1
  {fourier_path<PRECISION>(0.0, 1.0, 0.0, 0.0,
  fourier<PRECISION>(1, 1.0, c1)), 
  fourier_path<PRECISION>(0.0, 1.0, 1.0, 1.0,
  fourier<PRECISION>(1, 1.0, c2))}; // both vaild
  vector<fourier_path<PRECISION>> path2
  {fourier_path<PRECISION>(0.0, 1.0, 0.0, 0.0,
  fourier<PRECISION>(1, 4.0, c1)), 
  fourier_path<PRECISION>(0.0, 1.0, 1.0, 1.0,
  fourier<PRECISION>(1, 1.0, c2))}; // first path invalid
  vector<fourier_path<PRECISION>> path3
  {fourier_path<PRECISION>(0.0, 1.0, 0.0, 0.0,
  fourier<PRECISION>(1, 1.0, c1)), 
  fourier_path<PRECISION>(0.0, 1.0, 1.0, 1.0,
  fourier<PRECISION>(1, 2.0, c2))}; // second path invalid
  vector<fourier_path<PRECISION>> path4
  {fourier_path<PRECISION>(0.0, 1.0, 0.0, 0.0,
  fourier<PRECISION>(1, 4.0, c1)), 
  fourier_path<PRECISION>(0.0, 1.0, 1.0, 1.0,
  fourier<PRECISION>(1, 2.0, c2))}; // both path invalid
  action<PRECISION, fourier_path<PRECISION>, zero_lag<PRECISION>> tst;

  cout << "==========================================================" << endl;
  cout << "               Test action::is_vaild() routine            " << endl;
  cout << " Test 1. both element of path vector are vaild            " << endl;
  cout << " Test 2. first element of path vector is invaild          " << endl;
  cout << " Test 3. second element of path vector is invaild         " << endl;
  cout << " Test 4. both element of path vector is invaild           " << endl;
  cout << boolalpha;
  tst.update(path1);
  result = tst.is_vaild();
  sucess = sucess && (result == true);
  cout << "Test 1. vaildity of path: " << result << " it should be " << true << endl;
  tst.update(path2);
  result = tst.is_vaild();
  sucess = sucess && (result == false);
  cout << "Test 2. vaildity of path: " << result << " it should be " << false << endl;
  tst.update(path3);
  result = tst.is_vaild();
  sucess = sucess && (result == false);
  cout << "Test 3. vaildity of path: " << result << " it should be " << false << endl;
  tst.update(path4);
  result = tst.is_vaild();
  sucess = sucess && (result == false);
  cout << "Test 4. vaildity of path: " << result << " it should be " << false << endl;
  cout << "Test finished " << endl;
  if(sucess)
  return 0;
  cout << "Test failed" << endl;
  return -1;
}