.SUFFIXES: .cpp .o

CXX = g++
INC = -I src/include/
CXXFLAGS = -O3 $(INC)

OBJS = ./src/hw3.o ./src/support.o ./src/main.o
SRCS = ./src/hw3.cpp ./src/main.cpp ./src/support.cpp

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

./src/hw3.o: src/include/hw3.hpp
./src/main.o: src/include/hw3.hpp
./src/support.o: src/include/hw3.hpp
