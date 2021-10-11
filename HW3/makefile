.SUFFIXES: .cpp .o

CXX = g++
INC = -I src/include/
CXXFLAGS = -O2 -Wall -Wextra -pedantic $(INC) # default
# CXXFLAGS = -O3 $(INC) # optimization

OBJS = ./src/hw3.o ./src/main.o \
       ./src/kepler.o ./src/path.o ./src/fourier.o \
       ./src/move.o ./src/vector.o
SRCS = ./src/hw3.cpp ./src/main.cpp \
       ./src/kepler.cpp ./src/path.cpp ./src/fourier.cpp \
       ./src/move.cpp ./src/vector.cpp 

TARGET = ./bin/hw3

all : $(TARGET)

$(TARGET) : $(OBJS)
	$(CXX) -o $@ $(OBJS)

dep :
	makedepend $(INC) $(SRCS)

clean :
	rm -rf $(OBJS) $(TARGET)

new :
	$(MAKE) clean
	$(MAKE)
# DO NOT DELETE

./src/hw3.o: src/include/basic.hpp src/include/kepler.hpp
./src/main.o: src/include/hw3.hpp
./src/kepler.o: src/include/basic.hpp
./src/path.o: src/include/basic.hpp
./src/fourier.o: src/include/basic.hpp
