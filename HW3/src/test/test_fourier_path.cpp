/*!
 * @file test_fourier_path.cpp
 * @brief test fourier_path class routine 
 * @author pistack (Junho Lee)
 * @date 2021. 11. 5.
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
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
    PRECISION eps = 100*std::numeric_limits<PRECISION>::epsilon();
    bool sucess = true;
    vector<PRECISION> c {
        1.0, 0.0
    };
    vector<PRECISION> c2 {
        0.0, 1.0
    };
    PRECISION t0 = 0.0, t1 = 1.0;
    PRECISION p0 = 0.0, p1 = 0.0;
    PRECISION p0_2 = 0.0, p1_2 = 1.0;
    PRECISION add, scale;
    fourier<PRECISION> func1(1, 1.0, c);
    fourier<PRECISION> func2(1, 1.0, c2);
    fourier<PRECISION> func3(1, 2.0, c);
    fourier<PRECISION> func4(1, 2.0, c2);
    fourier_path<PRECISION> path(t0, t1, p0, p1);
    fourier_path<PRECISION> path2(t0, t1, p0_2, p1_2);

    cout << "==========================================================" << endl;
    cout << "Test fourier path class routine " << endl;
    cout << "Test 1. sin(2*pi*x) with boundary condition f(0)=f(1)=0 " << endl;
    cout << "Test 2. cos(2*pi*x) with boundary condition f(0)=f(1)=0 " << endl;
    cout << "Test 3. cos(pi*x) with boundary condition f(0)=0, f(1)=1 " << endl;
    cout << "Test 4. sin(pi*x) with boundary condition f(0)=0, f(1)=1 " << endl;
    cout << "Test 1." << endl;
    path.update(func1);
    add = path.get_adder();
    scale = path.get_scaler();
    sucess = sucess && (std::abs(add)<eps) && (std::abs(scale-1)<eps);
    cout << "Adder: " << path.get_adder() << " should be 0." << endl;
    cout << "Scaler: " << path.get_scaler() << " should be 1." << endl;
    cout << "Test 2." << endl;
    path.update(func2);
    add = path.get_adder();
    scale = path.get_scaler();
    sucess = sucess && (std::abs(add+1)<eps) && (std::abs(scale-1)<eps);
    cout << "Adder: " << path.get_adder() << " should be -1." << endl;
    cout << "Scaler: " << path.get_scaler() << " should be 1." << endl;
    cout << "Test 3." << endl;
    path2.update(func4);
    add = path2.get_adder();
    scale = path2.get_scaler();
    sucess = sucess && (std::abs(add-0.5)<eps) && (std::abs(scale+0.5)<eps);
    cout << "Adder: " << path2.get_adder() << " should be 0.5." << endl;
    cout << "Scaler: " << path2.get_scaler() << " should be -0.5." << endl;
    cout << "Test 4." << endl;
    path2.update(func3);
    add = path2.get_adder();
    scale = path2.get_scaler();
    sucess = sucess && (std::abs(add)<eps) && (std::abs(scale)<eps) && !path2.is_vaild();
    cout << boolalpha;
    cout << "Vaildity of path: " << path2.is_vaild() << " should be false. " << endl; 
    cout << "Adder: " << path2.get_adder() << " should be 0. " << endl;
    cout << "Scaler: " << path2.get_scaler() << " should be 0. " << endl;
    cout << "Test finished. " << endl;
    cout << "==========================================================" << endl;
    if(sucess)
    return 0;
    return -1;
}

