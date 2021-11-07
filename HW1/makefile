.SUFFIXES: .cpp .o

CXX = g++
INC = -I src/include/
#CPPFLAGS = -O3 -D PRECISION_LEVEL=0 $(INC) # for single precision
CPPFLAGS = -O3 -D PRECISION_LEVEL=1 $(INC) # for double precision

OBJS = ./src/main.o
SRCS = ./src/main.cpp

TARGET = ./bin/hw1

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

./src/main.o: src/include/hw1.hpp src/include/hw1.tpp
