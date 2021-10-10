.SUFFIXES: .cpp .o

CXX = g++
INC = -I src/include/
CXXFLAGS = -O3 $(INC)

OBJS = ./src/hw1.o ./src/main.o
SRCS = ./src/hw1.cpp ./src/main.cpp

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

./src/hw1.o: src/include/hw1.hpp
./src/main.o: src/include/hw1.hpp
