cmake_minimum_required(VERSION 3.0.0)
project(Compter_HW3 VERSION 0.1.0)
if(DEFINED PRECISION_LEVEL)
add_compile_definitions(PRECISION_LEVEL=${PRECISION_LEVEL})
else()
add_compile_definitions(PRECISION=1)
endif()

include(CTest)
enable_testing()
include_directories(src/include)

if(PRECISION_LEVEL EQUAL 0)
add_executable(hw3f
src/main.cpp src/include/action.hpp
src/include/fourier.hpp
src/include/fourier_path.hpp
src/include/mcm.hpp
src/include/action/action.tpp
src/include/fourier/fourier.tpp
src/include/fourier_path/fourier_path.tpp
src/include/mcm/mcm_basic.tpp
src/include/mcm/mcm_move.tpp
src/include/mcm/mcm_opt.tpp)
else()
add_executable(hw3d
src/main.cpp src/include/action.hpp
src/include/fourier.hpp
src/include/fourier_path.hpp
src/include/mcm.hpp
src/include/action/action.tpp
src/include/fourier/fourier.tpp
src/include/fourier_path/fourier_path.tpp
src/include/mcm/mcm_basic.tpp
src/include/mcm/mcm_move.tpp
src/include/mcm/mcm_opt.tpp)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
