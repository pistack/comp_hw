CXX = g++ # c++ complier name

# precision level 0: single, 1: double
# Monitor optimization process 0: no, 1: yes
# PATH_TYPE: PATH_TYPE_BEZIER: bezier curve
# PATH_TYPE: PATH_TYPE_FOURIER: fourier function

PRECISION_LEVEL = 1
MONITOR = 1 
PATH_TYPE = PATH_TYPE_BEZIER 

# C++ complier flag 
# example:
# DEBUG: -O -g -Wall -pedantic
# RELEASE: -O3

CXXFLAGS = -O3