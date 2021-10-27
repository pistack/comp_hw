.SUFFIXES: .cpp .o

CXX = g++
INC = -I src/include/
CXXFLAGS = -O -g -Wall -Werror -pedantic $(INC)
# CXXFLAGS = -O3 $(INC)
# CPPFLAGS = -D PRECISION=float # for float precision
CPPFLAGS = -D PRECISION=double # for double precision

OBJS = ./src/main.o
SRCS = ./src/main.cpp

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

./src/main.o: src/include/mcm.hpp src/include/fourier_path.hpp
./src/main.o: src/include/fourier.hpp src/include/fourier/fourier.tpp
./src/main.o: src/include/fourier_path/fourier_path.tpp
./src/main.o: src/include/action.hpp src/include/action/action.tpp
./src/main.o: src/include/mcm/mcm_basic.tpp src/include/mcm/mcm_move.tpp
./src/main.o: src/include/mcm/mcm_opt.tpp
