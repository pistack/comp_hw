#include <cmath>
#include <iostream>
#include <libpath/math_const.hpp>


int main(void)
{
    std::cout << libpath::PI<float>() << std::endl;
    std::cout << std::sin(libpath::PI<float>()) << std::endl;
    return 0;
}