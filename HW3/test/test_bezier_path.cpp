/*!
 * @file test_bezier_path.cpp
 * @brief test bezier_path class routine 
 * @author pistack (Junho Lee)
 * @date 2021. 11. 10.
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include "libpath/bezier_path.hpp"

#if PRECISION_LEVEL == 0
    #define PRECISION float
    #define DIGITS 6
#elif PRECISION_LEVEL == 1
    #define PRECISION double
    #define DIGITS 14
#endif

using namespace libpath;
using namespace std;

int main()
{
    bool sucess = true;
    PRECISION eps = 100*std::numeric_limits<PRECISION>::epsilon();
    vector<PRECISION> c1 = {1, 0}; // two control points -> order 1 bezier curve
    vector<PRECISION> c2 = {1, 2, 1}; // three control points -> order 2 bezier curve
    vector<PRECISION> c3 = {1, 2, 3, 1}; // four control points -> order 3 bezier curve
    vector<PRECISION> c4 = {17, 2, 3, 1}; // four control points -> order 3 bezier curve
    bezier_path<PRECISION> path(-1, 1, 0, 1);
    bezier_path<PRECISION> path2(-1, 1, 0.3, 0.5);
    bezier<PRECISION> fun1(1, c1); // f(t) = 1-t
    bezier<PRECISION> fun2(2, c2); // f(t) = -2t^2 + 2t+1 
    bezier<PRECISION> fun3(3, c3); // f(t) = -3t^3 + 3t + 1
    vector<PRECISION> t{
        -1, -1.0/3, 0, 1
    };
    vector<PRECISION> result;
    vector<vector<PRECISION>> result_eval
    {
        {0, 0, 0, 0},
        {0, 1, 5.0/4, 1},
        {1, 17.0/9,  17.0/8, 1}
    };
    vector<vector<PRECISION>> result_deriv
    {
        {0, 0, 0, 0},
        {2, 1, 0.5, -1},
        {3, 2, 3.0/4, -6},
    };  
    cout << "==================================================================" << endl;
    cout << "                   Test bezier_path class                         " << endl;
    cout << " initial condition of path:                                       " << endl;
    cout << " t_0: -1, t_1: 1, p_0 = 0, p_1 = 1                                " << endl;
    cout << " Test 1. order 1 bezier with control points: (1, 0)               " << endl;
    cout << " Test 2. order 2 bezier with control points: (1, 2, 1)            " << endl;
    cout << " Test 3. order 3 bezier with control points: (1, 2, 3, 1)         " << endl;
    cout << " initial condition of path2:                                      " << endl;
    // Test 1.
    cout << boolalpha;
    path.update(fun1); // invalid
    result = path.eval(t);
    for(int i=0; i<4; i++)
    {
        if(std::abs(result[i]-result_eval[0][i])>eps)
        sucess = false;
    }
    cout << "Test 1. order 1 bezier" << endl;
    cout << "Vaildity of bezier curve " << path.is_vaild() << endl;
    cout << "function value: " << endl;
    cout << "-1" << '\t' << "-1/3" << '\t' << "0" << '\t' << "1" << endl;
    for(int i=0; i<4; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = path.deriv(t);
    for(int i=0; i<4; i++)
    {
        if(std::abs(result[i]-result_deriv[0][i])>eps)
        sucess = false;
    }
    cout << "derivative: " << endl;
    cout << "-1" << '\t' << "-1/3" << '\t' << "0" << '\t' << "1" << endl;
    for(int i=0; i<4; ++i)
    cout << result[i] << '\t';
    cout << endl;
    // Test 2. 
    path.update(fun2); /// control point changes to (1, 2, 1) to (0, 2, 1)
    result = path.eval(t);
    for(int i=0; i<4; i++)
    {
        if(std::abs(result[i]-result_eval[1][i])>eps)
        sucess = false;
    }
    cout << "Test 2. order 2 bezier" << endl;
    cout << "Vaildity of bezier curve " << path.is_vaild() << endl;
    cout << "function value: " << endl;
    cout << "-1" << '\t' << "-1/3" << '\t' << "0" << '\t' << "1" << endl;
    for(int i=0; i<4; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = path.deriv(t);
    for(int i=0; i<4; i++)
    {
        if(std::abs(result[i]-result_deriv[1][i])>eps)
        sucess = false;
    }
    cout << "derivative: " << endl;
    cout << "-1" << '\t' << "-1/3" << '\t' << "0" << '\t' << "1" << endl;
    for(int i=0; i<4; ++i)
    cout << result[i] << '\t';
    cout << endl;
    // Test 3.
    result = fun3.eval(t);
    for(int i=0; i<4; i++)
    {
        if(std::abs(result[i]-result_eval[2][i])>eps)
        sucess = false;
    }
    cout << "Test 3. order 3 bezier" << endl;
    cout << "function value: " << endl;
    cout << "0" << '\t' << "1/3" << '\t' << "1/2" << '\t' << "1" << endl;
    for(int i=0; i<4; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = fun3.deriv(t);
    for(int i=0; i<4; i++)
    {
        if(std::abs(result[i]-result_deriv[2][i])>eps)
        sucess = false;
    }
    cout << "derivative: " << endl;
    cout << "0" << '\t' << "1/3" << '\t' << "1/2" << '\t' << "1" << endl;
    for(int i=0; i<4; ++i)
    cout << result[i] << '\t';
    cout << endl;
    cout << "Test finished" << endl;
    if(sucess)
    return 0;
    cout << "Test failed" << endl;
    return -1;
}