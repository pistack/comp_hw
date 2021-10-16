cmake_minimum_required(VERSION 3.0.0)
project(Compter_HW2 VERSION 0.1.0)

include(CTest)
enable_testing()
include_directories(src/include)
add_executable(hw2
src/main.cpp src/hw2.cpp 
src/include/hw2.hpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
