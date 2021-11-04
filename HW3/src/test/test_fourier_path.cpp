/*!
 * @file test_fourier_path.cpp
 * @brief test fourier_path class routine 
 * @author pistack (Junho Lee)
 * @date 2021. 11. 4.
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include "fourier_path.hpp"

using namespace std;
using namespace libfourier;

#if PRECISION_LEVEL == 0
    #define PRECISION float
    #define DIGITS 6
#elif PRECISION_LEVEL == 1
    #define PRECISION double
    #define DIGITS 14
#endif

int main(void)
{
    vector<PRECISION> c {
        0.0, 1.0
    };
    vector<PRECISION> c2 {
        1.0, 0.0
    };
    PRECISION t0 = 0.0, t1 = 1.0;
    PRECISION p0 = 0.0, p1 = 0.0;
    PRECISION p0_2 = 0.0, p1_2 = 1.0;
    fourier<PRECISION> func1(1, 1.0, c);
    fourier<PRECISION> func2(1, 1.0, c);
    fourier<PRECISION> func3(1, 2.0, c);
    fourier<PRECISION> func4(1, 2.0, c);

    cout << "Test fourier path class routine " << endl;
    cout << "Test 1. sin(2*pi*x) with boundary condition f(0)=f(1)=0 " << endl;
    cout << "Test 2. cos(2*pi*x) with boundary condition f(0)=f(1)=0 " << endl;
    cout << "Test 3. sin(pi*x) with boundary condition f(0)=0, f(1)=1 " << endl;
    cout << "Test 4. cos(pi*x) with boundary condition f(0)=0, f(1)=1 " << endl;
    return 0;
}

