/*!
 * @file test_fourier.cpp
 * @brief test fourier class routine 
 * @author pistack (Junho Lee)
 * @date 2021. 11. 6.
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include "fourier.hpp"

#if PRECISION_LEVEL == 0
    #define PRECISION float
    #define DIGITS 6
#elif PRECISION_LEVEL == 1
    #define PRECISION double
    #define DIGITS 14
#endif

using namespace libfourier;
using namespace std;

int main()
{
    bool sucess = true;
    PRECISION eps = 100*std::numeric_limits<PRECISION>::epsilon();
    vector<PRECISION> c1 = {1.0, 0.0};
    vector<PRECISION> c2 = {0.0, 1.0};
    vector<PRECISION> c3 = {1.0, 0.0, 1.0, 0.0};
    fourier<PRECISION> fun1(1, pi<PRECISION>, c1); // sin(2*x)
    fourier<PRECISION> fun2(1, pi<PRECISION>, c2); // cos(2*x)
    fourier<PRECISION> fun3(2, pi<PRECISION>, c3); // sin(2*x) + sin(4*x)
    vector<PRECISION> t{
        0.0, pi<PRECISION>/4, pi<PRECISION>/2, 3*pi<PRECISION>/4, pi<PRECISION>
    };
    vector<PRECISION> result;
    vector<vector<PRECISION>> result_eval
    {
        {0, 1, 0, -1, 0},
        {1, 0, -1, 0, 1},
        {0, 1, 0, -1, 0}
    };
    vector<vector<PRECISION>> result_deriv_sine
    {
        {2, 0, -2, 0, 2},
        {0, -4, 0, 4, 0},
        {-8, 0, 8, 0, -8},
        {0, 16, 0, -16, 0},
        {32, 0, -32, 0, 32}
    };  
    vector<vector<PRECISION>> result_deriv_cos
    {
        {0, -2, 0, 2, 0},
        {-4, 0, 4, 0, -4},
        {0, 8, 0, -8, 0},
        {16, 0, -16, 0, 16},
        {0, -32, 0, 32, 0}
    };

    cout << "==========================================================" << endl;
    cout << "               Test fourier class                         " << endl;
    cout << " Test 1. sin(2*x) at (0, pi/4, pi/2 3*pi/4, pi)           " << endl;
    cout << " Test 2. cos(2*x) at (0, pi/4, pi/2 3*pi/4, pi)           " << endl;
    cout << " Test 3. sin(2*x)+sin(4*x) at (0, pi/4, pi/2 3*pi/4, pi)  " << endl;
    cout << " Test 4. nth derivative of sin(2*x) for n=1,2,3,4,5       " << endl;
    cout << " Test 5. nth derivative of cos(2*x) for n=1,2,3,4,5       " << endl;
    result = fun1.eval(t);
    for(int i=0; i<5; i++)
    {
        if(std::abs(result[i]-result_eval[0][i])>eps)
        sucess = false;
    }
    cout << "Test 1. sin(2*x)" << endl;
    cout << 0 << '\t' << "pi/4" << '\t' << "pi/2" << '\t' << "3*pi/4" << '\t' << "pi" << endl;
    for(int i=0; i<5; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = fun2.eval(t);
    for(int i=0; i<5; i++)
    {
        if(std::abs(result[i]-result_eval[1][i])>eps)
        sucess = false;
    }
    cout << "Test 2. cos(2*x)" << endl;
    cout << 0 << '\t' << "pi/4" << '\t' << "pi/2" << '\t' << "3*pi/4" << '\t' << "pi" << endl;
    for(int i=0; i<5; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = fun3.eval(t);
    for(int i=0; i<5; i++)
    {
        if(std::abs(result[i]-result_eval[2][i])>eps)
        sucess = false;
    }
    cout << "Test 3. sin(2*x) + sin(4*x)" << endl;
    cout << 0 << '\t' << "pi/4" << '\t' << "pi/2" << '\t' << "3*pi/4" << '\t' << "pi" << endl;
    for(int i=0; i<5; ++i)
    cout << result[i] << '\t';
    cout << endl;
    cout << "Test 4. sin(2*x)" << endl;
    cout << "1st derivative of sin(2*x)" << endl;
    result = fun1.deriv(t);
    for(int i=0; i<5; i++)
    {
        if(std::abs(result[i]-result_deriv_sine[0][i])>eps)
        sucess = false;
    }
    cout << 0 << '\t' << "pi/4" << '\t' << "pi/2" << '\t' << "3*pi/4" << '\t' << "pi" << endl;
    for(int i=0; i<5; ++i)
    cout << result[i] << '\t';
    cout << endl;
    for(int i=2; i<6; ++i)
    {
        cout << i << "th derivative of sin(2*x)" << endl;
        result = fun1.nderiv(i, t);
        for(int j=0; j<5; ++j)
        {
            if(std::abs(result[j]-result_deriv_sine[i-1][j])>eps)
            sucess = false;
        }
        cout << 0 << '\t' << "pi/4" << '\t' << "pi/2" << '\t' << "3*pi/4" << '\t' << "pi" << endl;
        for(int j=0; j<5; ++j)
        cout << result[j] << '\t';
        cout << endl;
    }
    cout << "Test 5. cos(2*x)" << endl;
    cout << "1st derivative of cos(2*x)" << endl;
    result = fun2.deriv(t);
    for(int i=0; i<5; ++i)
    {
        if(std::abs(result[i]-result_deriv_cos[0][i])>eps)
        sucess = false;
    }
    cout << 0 << '\t' << "pi/4" << '\t' << "pi/2" << '\t' << "3*pi/4" << '\t' << "pi" << endl;
    for(int i=0; i<5; ++i)
    cout << result[i] << '\t';
    cout << endl;
    for(int i=2; i<6; ++i)
    {
        cout << i << "th derivative of cos(2*x)" << endl;
        result = fun2.nderiv(i, t);
        for(int j=0; j<5; ++j)
        {
            if(std::abs(result[j]-result_deriv_cos[i-1][j])>eps)
            sucess = false;
        }
        cout << 0 << '\t' << "pi/4" << '\t' << "pi/2" << '\t' << "3*pi/4" << '\t' << "pi" << endl;
        for(int i=0; i<5; ++i)
        cout << result[i] << '\t';
        cout << endl;
    }
    cout << "Test finished" << endl;
    if(sucess)
    return 0;
    return -1;
}