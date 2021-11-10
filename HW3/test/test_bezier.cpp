/*!
 * @file test_bezier.cpp
 * @brief test bezier class routine 
 * @author pistack (Junho Lee)
 * @date 2021. 11. 10.
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include "libpath/bezier.hpp"

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
    bezier<PRECISION> fun1(1, c1); // f(t) = 1-t
    bezier<PRECISION> fun2(2, c2); // f(t) = -2t^2 + 2t+1 
    bezier<PRECISION> fun3(3, c3); // f(t) = -3t^3 + 3t + 1
    vector<PRECISION> t{
        0, 1.0/3, 1.0/2, 1
    };
    vector<PRECISION> result;
    vector<vector<PRECISION>> result_eval
    {
        {1, 2.0/3, 1.0/2, 0},
        {1, 13.0/9, 3.0/2, 1},
        {1, 17.0/9,  17.0/8, 1}
    };
    vector<vector<PRECISION>> result_deriv
    {
        {-1, -1, -1, -1},
        {2, 2.0/3, 0, -2},
        {3, 2, 3.0/4, -6},
    };  
    cout << "==================================================================" << endl;
    cout << "                        Test bezier class                         " << endl;
    cout << " Test 1. order 1 bezier with control points: (1, 0)               " << endl;
    cout << " Test 2. order 2 bezier with control points: (1, 2, 1)            " << endl;
    cout << " Test 3. order 3 bezier with control points: (1, 2, 3, 1)         " << endl;
    // Test 1.
    result = fun1.eval(t);
    for(int i=0; i<4; i++)
    {
        if(std::abs(result[i]-result_eval[0][i])>eps)
        sucess = false;
    }
    cout << "Test 1. order 1 bezier" << endl;
    cout << "function value: " << endl;
    cout << "0" << '\t' << "1/3" << '\t' << "1/2" << '\t' << "1" << endl;
    for(int i=0; i<4; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = fun1.deriv(t);
    for(int i=0; i<4; i++)
    {
        if(std::abs(result[i]-result_deriv[0][i])>eps)
        sucess = false;
    }
    cout << "derivative: " << endl;
    cout << "0" << '\t' << "1/3" << '\t' << "1/2" << '\t' << "1" << endl;
    for(int i=0; i<4; ++i)
    cout << result[i] << '\t';
    cout << endl;
    // Test 2. 
    result = fun2.eval(t);
    for(int i=0; i<4; i++)
    {
        if(std::abs(result[i]-result_eval[1][i])>eps)
        sucess = false;
    }
    cout << "Test 2. order 2 bezier" << endl;
    cout << "function value: " << endl;
    cout << "0" << '\t' << "1/3" << '\t' << "1/2" << '\t' << "1" << endl;
    for(int i=0; i<4; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = fun2.deriv(t);
    for(int i=0; i<4; i++)
    {
        if(std::abs(result[i]-result_deriv[1][i])>eps)
        sucess = false;
    }
    cout << "derivative: " << endl;
    cout << "0" << '\t' << "1/3" << '\t' << "1/2" << '\t' << "1" << endl;
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