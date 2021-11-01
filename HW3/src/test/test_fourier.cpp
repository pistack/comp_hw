/*!
 * @file test_fourier.cpp
 * @brief test fourier class routine 
 * @author pistack (Junho Lee)
 * @date 2021. 11. 2.
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include "fourier.hpp"

#if PRECISION_LEVEL == 0
    #define PRECISION float
    #define DIGITS 6
#elif PRECISION_LEVEL == 1
    #define PRECISION double
    #define DIGITS 14
#endif

using namespace std;

const PRECISION pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062; // define pi

int main()
{
    vector<PRECISION> c1 = {1.0, 0.0};
    vector<PRECISION> c2 = {0.0, 1.0};
    fourier<PRECISION> fun1(1, pi, c1); // sin(2*x)
    fourier<PRECISION> fun2(1, pi, c2); // cos(2*x)
    vector<PRECISION> t{
        0.0, 0.25*pi, 0.5*pi, 0.75*pi, pi
    };
    vector<PRECISION> result(5,0);

    cout << "==========================================================" << endl;
    cout << "               Test fourier class                         " << endl;
    cout << " Test 1. sin(2*x) at (0, pi/4, pi/2 3*pi/4, pi)           " << endl;
    cout << " Test 2. cos(2*x) at (0, pi/4, pi/2 3*pi/4, pi)           " << endl;
    cout << " Test 3. nth derivative of sin(2*x) for n=1,2,3,4,5       " << endl;
    cout << " Test 4. nth derivative of cos(2*x) for n=1,2,3,4,5       " << endl;
    result = fun1.eval(t);
    cout << "Test 1. sin(2*x)" << endl;
    cout << 0 << '\t' << "pi/4" << '\t' << "pi/2" << '\t' << "3*pi/4" << '\t' << "pi" << endl;
    for(int i=0; i<5; ++i)
    cout << result[i] << '\t';
    cout << endl;
    result = fun2.eval(t);
    cout << "Test 2. cos(2*x)" << endl;
    cout << 0 << '\t' << "pi/4" << '\t' << "pi/2" << '\t' << "3*pi/4" << '\t' << "pi" << endl;
    for(int i=0; i<5; ++i)
    cout << result[i] << '\t';
    cout << endl;
    cout << "Test 3. sin(2*x)" << endl;
    cout << "1st derivative of sin(2*x)" << endl;
    result = fun1.deriv(t);
    cout << 0 << '\t' << "pi/4" << '\t' << "pi/2" << '\t' << "3*pi/4" << '\t' << "pi" << endl;
    for(int i=0; i<5; ++i)
    cout << result[i] << '\t';
    cout << endl;
    for(int i=2; i<6; ++i)
    {
        cout << i << "th derivative of sin(2*x)" << endl;
        result = fun1.nderiv(i, t);
        cout << 0 << '\t' << "pi/4" << '\t' << "pi/2" << '\t' << "3*pi/4" << '\t' << "pi" << endl;
        for(int i=0; i<5; ++i)
        cout << result[i] << '\t';
        cout << endl;
    }
    cout << "Test 4. cos(2*x)" << endl;
    cout << "1st derivative of cos(2*x)" << endl;
    result = fun2.deriv(t);
    cout << 0 << '\t' << "pi/4" << '\t' << "pi/2" << '\t' << "3*pi/4" << '\t' << "pi" << endl;
    for(int i=0; i<5; ++i)
    cout << result[i] << '\t';
    cout << endl;
    for(int i=2; i<6; ++i)
    {
        cout << i << "th derivative of cos(2*x)" << endl;
        result = fun2.nderiv(i, t);
        cout << 0 << '\t' << "pi/4" << '\t' << "pi/2" << '\t' << "3*pi/4" << '\t' << "pi" << endl;
        for(int i=0; i<5; ++i)
        cout << result[i] << '\t';
        cout << endl;
    }
    cout << "Test finished" << endl;
    return 0;
}