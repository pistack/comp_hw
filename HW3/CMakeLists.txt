cmake_minimum_required(VERSION 3.0.0)
project(Compter_HW3 VERSION 0.1.0)

include(CTest)
enable_testing()
include_directories(src/include)
add_library(fourier STATIC
src/fourier.cpp src/fourier_path.cpp
src/action.cpp src/include/fourier.hpp
src/include/fourier_path.hpp
src/include/action.hpp)
add_executable(hw3
src/main.cpp src/hw3.cpp src/move.cpp
src/include/hw3.hpp)
target_link_libraries(hw3 fourier)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
