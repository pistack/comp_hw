cmake_minimum_required(VERSION 3.0.0)
project(Compter_HW3 VERSION 0.1.0)
add_compile_definitions(PRECISION_LEVEL=1)
add_compile_definitions(MONITOR=0)
if(DEFINED PRECISION_LEVEL)
add_compile_definitions(PRECISION_LEVEL=${PRECISION_LEVEL})
endif()

if(DEFINED MONITOR)
if(MONITOR EQUAL "yes" OR MONITOR EQUAL 1)
add_compile_definitions(MONITOR=1)
endif()
endif()

include_directories(include)

if(PRECISION_LEVEL EQUAL 0)
add_executable(hw3f
src/main.cpp include/libpath/action.hpp
include/libpath/fourier.hpp
include/libpath/fourier_path.hpp
include/libpath/bezier.hpp
include/libpath/bezier_path.hpp
include/libmcm/mcm.hpp
include/libpath/node_weight_table.hpp
include/libpath/action/action.tpp
include/libpath/fourier/fourier.tpp
include/libpath/fourier_path/fourier_path.tpp
include/libpath/bezier/bezier.hpp
include/libpath/bezier_path/bezier_path.hpp
include/libmcm/mcm/mcm_basic.tpp
include/libmcm/mcm/mcm_move.tpp
include/libmcm/mcm/mcm_opt.tpp)
else()
add_executable(hw3d
src/main.cpp include/libpath/action.hpp
include/libpath/fourier.hpp
include/libpath/fourier_path.hpp
include/libpath/bezier.hpp
include/libpath/bezier_path.hpp
include/libmcm/mcm.hpp
include/libpath/node_weight_table.hpp
include/libpath/action/action.tpp
include/libpath/fourier/fourier.tpp
include/libpath/fourier_path/fourier_path.tpp
include/libpath/bezier/bezier.hpp
include/libpath/bezier_path/bezier_path.hpp
include/libmcm/mcm/mcm_basic.tpp
include/libmcm/mcm/mcm_move.tpp
include/libmcm/mcm/mcm_opt.tpp)
endif()

include(CTest)
enable_testing()
add_test()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
