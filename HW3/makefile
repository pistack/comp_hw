.SUFFIXES: .cpp .o

CXX = g++
AR = ar
INC = -I src/include/
CXXFLAGS = -O -g -Wall -Wextra -pedantic $(INC) #debug
#CXXFLAGS = -O2 -Wall -Wextra -pedantic $(INC) # default
#CXXFLAGS = -O3 $(INC) # optimization

OBJS = ./src/hw3.o ./src/main.o ./src/move.o
SRCS = ./src/hw3.cpp ./src/main.cpp \
       ./src/action.cpp ./src/fourier_path.cpp ./src/fourier.cpp \
       ./src/move.cpp

TARGET = hw3 lib

all : $(TARGET)

.PHONY: all

lib : ./src/fourier_path.o ./src/fourier.o ./src/action.o
	$(AR) rvs ./src/libfourier.a ./src/fourier_path.o ./src/fourier.o ./src/action.o

hw3 : lib $(OBJS)
	$(CXX) -o ./bin/hw3  $(OBJS) ./src/libfourier.a

dep :
	makedepend $(INC) $(SRCS)

clean :
	rm -rf $(OBJS) ./src/fourier_path.o ./src/fourier.o ./src/action.o

distclean :
	rm -rf $(OBJS) ./src/libfourier.a ./bin/hw3 ./src/fourier_path.o ./src/fourier.o ./src/action.o

new :
	$(MAKE) distclean
	$(MAKE)
# DO NOT DELETE

./src/hw3.o: src/include/fourier_path.hpp src/include/fourier.hpp
./src/hw3.o: src/include/action.hpp src/include/hw3.hpp
./src/main.o: src/include/hw3.hpp
./src/action.o: src/include/action.hpp src/include/fourier_path.hpp
./src/action.o: src/include/fourier.hpp
./src/fourier_path.o: src/include/fourier_path.hpp src/include/fourier.hpp
./src/fourier.o: src/include/fourier.hpp
./src/move.o: src/include/hw3.hpp
