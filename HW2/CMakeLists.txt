cmake_minimum_required(VERSION 3.0.0)
project(Compter_HW2 VERSION 0.1.0)
if(DEFINED PRECISION_LEVEL)
add_compile_definitions(PRECISION_LEVEL=${PRECISION_LEVEL})
else()
add_compile_definitions(PRECISION_LEVEL=1)
endif()

include(CTest)
enable_testing()
include_directories(src/include)

if(PRECISION_LEVEL EQUAL 0)
add_executable(hw2f
src/main.cpp src/include/hw2.tpp 
src/include/hw2.hpp)
else()
add_executable(hw2d
src/main.cpp src/include/hw2.tpp 
src/include/hw2.hpp)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
