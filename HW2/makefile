.SUFFIXES: .cpp .o

CXX = g++
INC = -I src/include/
CXXFLAGS = -O3 $(INC)
# CPPFLAGS = -D PRECISION_LEVEL=0 # for float precision
CPPFLAGS = -D PRECISION_LEVEL=1 # for double precision

OBJS = ./src/main.o
SRCS = ./src/main.cpp

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

./src/main.o: src/include/hw2.hpp src/include/hw2.tpp
