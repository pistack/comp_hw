.SUFFIXES: .cpp .o

CXX = g++
AR = ar
INC = -I src/include/
#CXXFLAGS = -O -g -Wall -Wextra -pedantic $(INC) #debug
#CXXFLAGS = -O2 -Wall -Wextra -pedantic $(INC) # default
CXXFLAGS = -O3 $(INC) # optimization
#CXXFLAGS = -O3 -mavx2 $(INC) # optimization for avx2 architecture

OBJS = ./src/main.o \
       ./src/mcm_basic.o ./src/mcm_move.o ./src/mcm_opt.o \
	   ./src/action.o ./src/fourier_path.o ./src/fourier.o
SRCS = ./src/main.cpp \
       ./src/mcm_basic.cpp ./src/mcm_move.cpp ./src/mcm_opt.cpp \
	   ./src/action.cpp ./src/fourier_path.cpp ./src/fourier.cpp

TARGET = hw3 lib

all : $(TARGET)

.PHONY: all

lib : ./src/fourier_path.o ./src/fourier.o ./src/action.o \
      ./src/mcm_basic.o ./src/mcm_opt.o ./src/mcm_move.o
	$(AR) rvs ./src/libfourier.a ./src/fourier_path.o ./src/fourier.o ./src/action.o
	$(AR) rvs ./src/libmcm.a ./src/mcm_basic.o ./src/mcm_opt.o ./src/mcm_move.o

hw3 : lib ./src/main.o
	$(CXX) -o ./bin/hw3  ./src/main.o ./src/libmcm.a ./src/libfourier.a

dep :
	makedepend $(INC) $(SRCS)

clean :
	rm -rf $(OBJS) 

distclean :
	rm -rf $(OBJS) ./src/libmcm.a ./src/libfourier.a ./bin/hw3

new :
	$(MAKE) distclean
	$(MAKE)
# DO NOT DELETE

./src/main.o: src/include/mcm.hpp src/include/fourier_path.hpp
./src/main.o: src/include/fourier.hpp src/include/action.hpp
./src/mcm_basic.o: src/include/mcm.hpp src/include/fourier_path.hpp
./src/mcm_basic.o: src/include/fourier.hpp src/include/action.hpp
./src/mcm_move.o: src/include/mcm.hpp src/include/fourier_path.hpp
./src/mcm_move.o: src/include/fourier.hpp src/include/action.hpp
./src/mcm_opt.o: src/include/mcm.hpp src/include/fourier_path.hpp
./src/mcm_opt.o: src/include/fourier.hpp src/include/action.hpp
./src/action.o: src/include/action.hpp src/include/fourier_path.hpp
./src/action.o: src/include/fourier.hpp
./src/fourier_path.o: src/include/fourier_path.hpp src/include/fourier.hpp
./src/fourier.o: src/include/fourier.hpp
