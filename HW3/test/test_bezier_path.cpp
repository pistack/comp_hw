/*!
 * @file test_bezier_path.cpp
 * @brief test bezier_path class routine 
 * @author pistack (Junho Lee)
 * @date 2021. 11. 12.
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include "libpath/bezier_path.hpp"
#include "test.hpp"

using namespace libpath;
using namespace std;

namespace test {
int test_bezier_path()
{
    bool sucess = true;
    PRECISION eps = 100*std::numeric_limits<PRECISION>::epsilon();
    int test_num = 0; // test number
    vector<PRECISION> c1 = {1, 0}; // two control points -> order 1 bezier curve
    vector<PRECISION> c2 = {1, 2}; 
    vector<PRECISION> c3 = {0, 1}; 
    bezier_path<PRECISION> path(-1, 1, 0, 1);
    bezier_path<PRECISION> path2(-1, 1, 2, 3);
    bezier_path<PRECISION> path3(-1, 1, 1, 0);
    bezier<PRECISION> fun1(1, c1); // f(t) = 1-t
    bezier<PRECISION> fun2(1, c2); // f(t) = -2t^2 + 2t+1 
    bezier<PRECISION> fun3(1, c3); // f(t) = -3t^3 + 3t + 1
    vector<PRECISION> t{
        -1, 1
    };
    vector<PRECISION> result;
    vector<vector<PRECISION>> result_eval
    {
        {0, 0}, {0, 1}, {0, 1},
        {0, 0}, {2, 3}, {0, 0},
        {1, 0}, {1, 0}, {0, 0}
    };
    vector<vector<PRECISION>> result_deriv
    {
        {0, 0}, {1./2, 1./2}, {1./2, 1./2},
        {0, 0}, {1./2, 1./2}, {0, 0},
        {-1./2, -1./2}, {-1./2, -1./2}, {0, 0}
    };  
    cout << "==================================================================" << endl;
    cout << "                   Test bezier_path class                         " << endl;
    cout << " initial condition of path:                                       " << endl;
    cout << " t_0: -1, t_1: 1, p_0 = 0, p_1 = 1                                " << endl;
    cout << " Test 1. order 1 bezier with control points: (1, 0)               " << endl;
    cout << " Test 2. order 1 bezier with control points: (1, 2)               " << endl;
    cout << " Test 3. order 1 bezier with control points: (0, 1)               " << endl;
    cout << " initial condition of path2:                                      " << endl;
    cout << " t_0: -1, t_1: 1, p_0 = 2, p_1 = 3                                " << endl;
    cout << " Test 4. order 1 bezier with control points: (1, 0)               " << endl;
    cout << " Test 5. order 1 bezier with control points: (1, 2)               " << endl;
    cout << " Test 6. order 1 bezier with control points: (0, 1)               " << endl;
    cout << " initial condition of path3:                                      " << endl;
    cout << " t_0: -1, t_1: 1, p_0 = 1, p_1 = 0                                " << endl;
    cout << " Test 7. order 1 bezier with control points: (1, 0)               " << endl;
    cout << " Test 8. order 1 bezier with control points: (1, 2)               " << endl;
    cout << " Test 9. order 1 bezier with control points: (0, 1)               " << endl;
    cout << "========================test start!===============================" << endl;
    cout << boolalpha;
    // Test 1.
    path.update(fun1); // invalid
    result = path.eval(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_eval[test_num][i])>eps)
        sucess = false;
    }
    cout << "Test 1" << endl;
    cout << "Vaildity of bezier curve " << path.is_vaild() << endl;
    cout << "function value: " << endl;
    cout << "-1" << '\t' <<  "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = path.deriv(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_deriv[test_num][i])>eps)
        sucess = false;
    }
    cout << "derivative: " << endl;
    cout << "-1" << '\t' << "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    ++test_num;
    // Test 2.
    path.update(fun2); // invalid
    result = path.eval(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_eval[test_num][i])>eps)
        sucess = false;
    }
    cout << "Test 2" << endl;
    cout << "Vaildity of bezier curve " << path.is_vaild() << endl;
    cout << "function value: " << endl;
    cout << "-1" << '\t' <<  "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = path.deriv(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_deriv[test_num][i])>eps)
        sucess = false;
    }
    cout << "derivative: " << endl;
    cout << "-1" << '\t' << "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    ++test_num;
    // Test 3.
    path.update(fun3); // invalid
    result = path.eval(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_eval[test_num][i])>eps)
        sucess = false;
    }
    cout << "Test 3" << endl;
    cout << "Vaildity of bezier curve " << path.is_vaild() << endl;
    cout << "function value: " << endl;
    cout << "-1" << '\t' <<  "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = path.deriv(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_deriv[test_num][i])>eps)
        sucess = false;
    }
    cout << "derivative: " << endl;
    cout << "-1" << '\t' << "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    ++test_num;
    // Test 4.
    path2.update(fun1); // invalid
    result = path2.eval(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_eval[test_num][i])>eps)
        sucess = false;
    }
    cout << "Test 4" << endl;
    cout << "Vaildity of bezier curve " << path2.is_vaild() << endl;
    cout << "function value: " << endl;
    cout << "-1" << '\t' <<  "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = path2.deriv(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_deriv[test_num][i])>eps)
        sucess = false;
    }
    cout << "derivative: " << endl;
    cout << "-1" << '\t' << "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    ++test_num;
    // Test 5.
    path2.update(fun2); // invalid
    result = path2.eval(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_eval[test_num][i])>eps)
        sucess = false;
    }
    cout << "Test 5" << endl;
    cout << "Vaildity of bezier curve " << path2.is_vaild() << endl;
    cout << "function value: " << endl;
    cout << "-1" << '\t' <<  "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = path2.deriv(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_deriv[test_num][i])>eps)
        sucess = false;
    }
    cout << "derivative: " << endl;
    cout << "-1" << '\t' << "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    ++test_num;
    // Test 6.
    path2.update(fun3); // invalid
    result = path2.eval(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_eval[test_num][i])>eps)
        sucess = false;
    }
    cout << "Test 6" << endl;
    cout << "Vaildity of bezier curve " << path2.is_vaild() << endl;
    cout << "function value: " << endl;
    cout << "-1" << '\t' <<  "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = path2.deriv(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_deriv[test_num][i])>eps)
        sucess = false;
    }
    cout << "derivative: " << endl;
    cout << "-1" << '\t' << "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    ++test_num;
    // Test 7.
    path3.update(fun1); // invalid
    result = path3.eval(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_eval[test_num][i])>eps)
        sucess = false;
    }
    cout << "Test 7" << endl;
    cout << "Vaildity of bezier curve " << path3.is_vaild() << endl;
    cout << "function value: " << endl;
    cout << "-1" << '\t' <<  "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = path3.deriv(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_deriv[test_num][i])>eps)
        sucess = false;
    }
    cout << "derivative: " << endl;
    cout << "-1" << '\t' << "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    ++test_num;
    // Test 8.
    path3.update(fun2); // invalid
    result = path3.eval(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_eval[test_num][i])>eps)
        sucess = false;
    }
    cout << "Test 8" << endl;
    cout << "Vaildity of bezier curve " << path3.is_vaild() << endl;
    cout << "function value: " << endl;
    cout << "-1" << '\t' <<  "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = path3.deriv(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_deriv[test_num][i])>eps)
        sucess = false;
    }
    cout << "derivative: " << endl;
    cout << "-1" << '\t' << "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    ++test_num;
    // Test 9.
    path3.update(fun3); // invalid
    result = path3.eval(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_eval[test_num][i])>eps)
        sucess = false;
    }
    cout << "Test 9" << endl;
    cout << "Vaildity of bezier curve " << path3.is_vaild() << endl;
    cout << "function value: " << endl;
    cout << "-1" << '\t' <<  "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = path3.deriv(t);
    for(int i=0; i<2; i++)
    {
        if(std::abs(result[i]-result_deriv[test_num][i])>eps)
        sucess = false;
    }
    cout << "derivative: " << endl;
    cout << "-1" << '\t' << "1" << endl;
    for(int i=0; i<2; ++i)
    cout << result[i] << '\t';
    cout << endl;
    ++test_num;

    cout << "Test finished" << endl;
    if(sucess)
    return 0;
    cout << "Test failed" << endl;
    return -1;
}
}