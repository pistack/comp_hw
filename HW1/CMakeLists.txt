cmake_minimum_required(VERSION 3.0.0)
project(Compter_HW1 VERSION 0.1.0)
if(DEFINED PRECISION)
add_compile_definitions(PRECISION=${PRECISION})
else()
add_compile_definitions(PRECISION=double)
endif()
include(CTest)
enable_testing()
include_directories(src/include)
if(PRECISION EQUAL "float")
add_executable(hw1f
src/main.cpp src/include/hw1.tpp 
src/include/hw1.hpp)
else()
add_executable(hw1d
src/main.cpp src/include/hw1.tpp 
src/include/hw1.hpp)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
