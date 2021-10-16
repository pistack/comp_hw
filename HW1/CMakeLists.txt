cmake_minimum_required(VERSION 3.0.0)
project(Compter_HW1 VERSION 0.1.0)

include(CTest)
enable_testing()
include_directories(src/include)
add_executable(hw1
src/main.cpp src/hw1.cpp 
src/include/hw1.hpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
