/*!
 * @file test.cpp
 * @brief test libpath
 * @author pistack (Junho Lee)
 * @date 2021. 11. 12.
 */

#include <chrono>
#include <vector>
#include <iostream>
#include <iomanip>
#include "libpath/math_const.hpp"
#include "test.hpp"

int main(void)
{
    // test report
    std::vector<bool> tst_result(7, false);
    std::vector<double> tst_elapsed(7, 0);

    //elapsed time
    std::chrono::steady_clock::time_point start;
    std::chrono::steady_clock::time_point end;

    // tst number
    int tst_num = 0;

    std::cout << libpath::PI<PRECISION>() << std::endl;
    std::cout << std::boolalpha;
    std::cout << "Testing libpath " << std::endl;
    std::cout << "Test 1. test bezier class " << std::endl;
    start = std::chrono::steady_clock::now();
    tst_result[tst_num] = !bool(test::test_bezier());
    end = std::chrono::steady_clock::now();
    tst_elapsed[tst_num] = \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000000.0;
    ++tst_num;
    std::cout << "Test 2. test bezier_path class " << std::endl;
    start = std::chrono::steady_clock::now();
    tst_result[tst_num] = !bool(test::test_bezier_path());
    end = std::chrono::steady_clock::now();
    tst_elapsed[tst_num] = \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000000.0;
    ++tst_num;
    std::cout << "Test 3. test fourier class " << std::endl;
    start = std::chrono::steady_clock::now();
    tst_result[tst_num] = !bool(test::test_fourier());
    end = std::chrono::steady_clock::now();
    tst_elapsed[tst_num] = \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000000.0;
    ++tst_num;
    std::cout << "Test 4. test fourier_path class " << std::endl;
    start = std::chrono::steady_clock::now();
    tst_result[tst_num] = !bool(test::test_fourier_path());
    end = std::chrono::steady_clock::now();
    tst_elapsed[tst_num] = \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000000.0;
    ++tst_num;
    std::cout << "Test action class " << std::endl;
    std::cout << "Test 5. test is_vaild() routine " << std::endl;
    start = std::chrono::steady_clock::now();
    tst_result[tst_num] = !bool(test::test_action_vaildity());
    end = std::chrono::steady_clock::now();
    tst_elapsed[tst_num] = \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000000.0;
    ++tst_num;
    std::cout << "Test 6. test eval() routine with various simple lagrangian " << std::endl;
    start = std::chrono::steady_clock::now();
    tst_result[tst_num] = !bool(test::test_action_simple());
    end = std::chrono::steady_clock::now();
    tst_elapsed[tst_num] = \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000000.0;
    ++tst_num;
    std::cout << "Test 7. test eval() routine with kepler lagrangian " << std::endl;
    start = std::chrono::steady_clock::now();
    tst_result[tst_num] = !bool(test::test_action_kepler());
    end = std::chrono::steady_clock::now();
    tst_elapsed[tst_num] = \
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000000.0;
    ++tst_num;
    tst_num = 0;
    // report test result
    std::cout << "Test Report " << std::endl;
    std::cout << "Test " << tst_num+1 << ". state: " << tst_result[tst_num] << " elapsed time: " << \
    tst_elapsed[tst_num] << " s" << std::endl;
    ++tst_num;
    std::cout << "Test " << tst_num+1 << ". state: " << tst_result[tst_num] << " elapsed time: " << \
    tst_elapsed[tst_num] << " s" << std::endl;
    ++tst_num;
    std::cout << "Test " << tst_num+1 << ". state: " << tst_result[tst_num] << " elapsed time: " << \
    tst_elapsed[tst_num] << " s" << std::endl;
    ++tst_num;
    std::cout << "Test " << tst_num+1 << ". state: " << tst_result[tst_num] << " elapsed time: " << \
    tst_elapsed[tst_num] << " s" << std::endl;
    ++tst_num;
    std::cout << "Test " << tst_num+1 << ". state: " << tst_result[tst_num] << " elapsed time: " << \
    tst_elapsed[tst_num] << " s" << std::endl;
    ++tst_num;
    std::cout << "Test " << tst_num+1 << ". state: " << tst_result[tst_num] << " elapsed time: " << \
    tst_elapsed[tst_num] << " s" << std::endl;
    ++tst_num;
    std::cout << "Test " << tst_num+1 << ". state: " << tst_result[tst_num] << " elapsed time: " << \
    tst_elapsed[tst_num] << " s" << std::endl;
    ++tst_num;
    std::cout << "Test finished! " << std::endl;
    return 0;
}