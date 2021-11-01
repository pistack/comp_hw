/*!
 * @file test_action_vaildity.cpp
 * @brief test action::is_vaild() routine 
 * @author pistack (Junho Lee)
 * @date 2021. 11. 2.
 */

#include <chrono>
#include <cmath>
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

int main()
{
    vector<PRECISION> c1 = {1.0, 0.0};
    vector<PRECISION> c2 = {0.0, 1.0};
    
    return 0;
}