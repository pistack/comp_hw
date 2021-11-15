cmake_minimum_required(VERSION 3.0.0)
project(Compter_HW3 VERSION 0.1.0)

if(DEFINED PRECISION_LEVEL)
add_compile_definitions(PRECISION_LEVEL=${PRECISION_LEVEL})
endif()

if(DEFINED MONITOR)
if(MONITOR STREQUAL "yes" OR MONITOR EQUAL 1)
add_compile_definitions(MONITOR=1)
endif()
endif()

if(DEFINED PATH_TYPE)
if((PATH_TYPE STREQUAL "bezier"))
add_compile_definitions(PATH_TYPE_BEZIER=1)
endif()
if(PATH_TYPE STREQUAL "fourier")
add_compile_definitions(PATH_TYPE_FOURIER=1)
endif()
endif()

if(NOT DEFINED PRECISION_LEVEL)
set(PRECISION_LEVEL 1)
endif()
if(NOT DEFINED MONITOR)
set(MONITOR 0)
endif()
if(NOT DEFINED PATH_TYPE)
set(PATH_TYPE "fourier")
endif()

include_directories(include)

# default include source
set(default_include
"include/libpath/math_const.hpp"
"include/libpath/action.hpp"
"include/libmcm/mcm.hpp"
"include/libpath/node_weight_table.hpp"
"include/libpath/math_const/math_const.tpp"
"include/libpath/action/action.tpp"
"include/libmcm/mcm/mcm_basic.tpp"
"include/libmcm/mcm/mcm_move.tpp"
"include/libmcm/mcm/mcm_opt.tpp"
)

# include source for beizer curve
set(beizer_include 
"include/libpath/bezier.hpp"
"include/libpath/bezier_path.hpp"
"include/libmcm/mcm_bezier.hpp"
"include/libpath/bezier/bezier.tpp"
"include/libpath/bezier_path/bezier_path.tpp"
"include/libmcm/mcm_bezier/mcm_bezier.tpp"
)

# include source for fourier function
set(fourier_include 
"include/libpath/fourier.hpp"
"include/libpath/fourier_path.hpp"
"include/libmcm/mcm_fourier.hpp"
"include/libpath/fourier/fourier.tpp"
"include/libpath/fourier_path/fourier_path.tpp"
"include/libmcm/mcm_fourier/mcm_fourier.tpp"
)

# fourier
if((PRECISION_LEVEL EQUAL 0) AND (PATH_TYPE STREQUAL "fourier") )
add_executable(hw3f_f
src/main.cpp ${default_include} ${fourier_include})
endif()
if((PRECISION_LEVEL EQUAL 1) AND (PATH_TYPE STREQUAL "fourier") )
add_executable(hw3f_d
src/main.cpp ${default_include} ${fourier_include})
endif()

# bezier
if((PRECISION_LEVEL EQUAL 0) AND (PATH_TYPE STREQUAL "bezier") )
add_executable(hw3b_f
src/main.cpp ${default_include} ${bezier_include})
endif()
if((PRECISION_LEVEL EQUAL 1) AND (PATH_TYPE STREQUAL "bezier") )
add_executable(hw3b_d
src/main.cpp ${default_include} ${bezier_include})
endif()

include(CTest)
enable_testing()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
