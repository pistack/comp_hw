.SUFFIXES: .cpp .o

CXX = g++
INC = -I src/include/
CXXFLAGS = -O3 $(INC)
# CPPFLAGS = -D PRECISION=float # for float precision
CPPFLAGS = -D PRECISION=double # for double precision

OBJS = ./src/hw2.o ./src/main.o
SRCS = ./src/hw2.cpp ./src/main.cpp

TARGET = ./bin/hw2

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

./src/hw2.o: src/include/hw2.hpp
./src/main.o: src/include/hw2.hpp
