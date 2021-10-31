.SUFFIXES: .cpp .o

CXX = g++
INC = -I src/include/
# CXXFLAGS = -O -g -Wall -Werror -pedantic $(INC)
# CPPFLAGS = -O3 -D PRECISION_LEVEL=0 $(INC) # for single precision
CPPFLAGS = -O3 -D PRECISION_LEVEL=1 $(INC) # for double precision

OBJS = ./src/main.o 
TEST_OBJS = ./src/test/test_action_kepler.o ./src/test/test_action_simple.o
SRCS = ./src/main.cpp ./src/test/test_action_kepler.cpp \
./src/test/test_action_simple.cpp

TARGET = ./bin/hw3
TEST = test_action_kepler test_action_simple

all : $(TARGET)

$(TARGET) : $(OBJS)
	$(CXX) -o $@ $(OBJS)

dep :
	makedepend $(INC) $(SRCS)

clean :
	rm -rf $(OBJS) $(TARGET) $(TEST) $(TEST_OBJS)

test : $(TEST)
	./test_action_kepler
	./test_action_simple
	rm -rf $(TEST) $(TEST_OBJS)

test_action_kepler: ./src/test/test_action_kepler.o
	$(CXX) -o $@ ./src/test/test_action_kepler.o

test_action_simple: ./src/test/test_action_simple.o
	$(CXX) -o $@ ./src/test/test_action_simple.o

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
./src/test/test_action_kepler.o: src/include/action.hpp
./src/test/test_action_kepler.o: src/include/fourier_path.hpp
./src/test/test_action_kepler.o: src/include/fourier.hpp
./src/test/test_action_kepler.o: src/include/fourier/fourier.tpp
./src/test/test_action_kepler.o: src/include/fourier_path/fourier_path.tpp
./src/test/test_action_kepler.o: src/include/action/action.tpp
./src/test/test_action_simple.o: src/include/action.hpp
./src/test/test_action_simple.o: src/include/fourier_path.hpp
./src/test/test_action_simple.o: src/include/fourier.hpp
./src/test/test_action_simple.o: src/include/fourier/fourier.tpp
./src/test/test_action_simple.o: src/include/fourier_path/fourier_path.tpp
./src/test/test_action_simple.o: src/include/action/action.tpp
