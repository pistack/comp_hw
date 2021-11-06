.SUFFIXES: .cpp .o

CXX = g++
INC = -I src/include/
# CPPFLAGS = -O -g -Wall -Werror -pedantic -D MONITOR=0 -D PRECISION_LEVEL=1 $(INC)
CPPFLAGS = -O3 -D MONITOR=1 -D PRECISION_LEVEL=0 $(INC) # for single precision
# CPPFLAGS = -O3 -D MONITOR=1 -D PRECISION_LEVEL=1 $(INC) # for double precision
# CPPFLAGS = -O3 -D MONITOR=0 -D PRECISION_LEVEL=0 $(INC) # for single precision
# CPPFLAGS = -O3 -D MONITOR=0 -D PRECISION_LEVEL=1 $(INC) # for double precision

OBJS = ./src/main.o 
TEST1_OBJS = ./src/test/test_action_vaildity.o
TEST2_OBJS = ./src/test/test_action_simple.o
TEST3_OBJS = ./src/test/test_action_kepler.o 
TEST4_OBJS = ./src/test/test_fourier.o 
TEST5_OBJS = ./src/test/test_fourier_path.o
TEST_OBJS = $(TEST1_OBJS) $(TEST2_OBJS) $(TEST3_OBJS) $(TEST4_OBJS) $(TEST5_OBJS)
SRCS = ./src/main.cpp ./src/test/test_action_kepler.cpp \
./src/test/test_action_simple.cpp ./src/test/test_action_vaildity.cpp \
./src/test/test_fourier.cpp ./src/test/test_fourier_path.cpp

TARGET = ./bin/hw3
TEST1 = test_action_vaildity
TEST2 = test_action_simple
TEST3 = test_action_kepler
TEST4 = test_fourier
TEST5 = test_fourier_path
TEST = $(TEST1) $(TEST2) $(TEST3) $(TEST4) $(TEST5)


all : $(TARGET)

$(TARGET) : $(OBJS)
	$(CXX) -o $@ $(OBJS)

dep :
	makedepend $(INC) $(SRCS)

clean :
	rm -rf $(OBJS) $(TARGET) $(TEST) $(TEST_OBJS)

test : $(TEST)
	./$(TEST1)
	./$(TEST2)
	./$(TEST3)
	./$(TEST4)
	./$(TEST5)
	rm -rf $(TEST) $(TEST_OBJS)

$(TEST1) : $(TEST1_OBJS)
	$(CXX) -o $@ $(TEST1_OBJS)

$(TEST2) : $(TEST2_OBJS)
	$(CXX) -o $@ $(TEST2_OBJS)

$(TEST3) : $(TEST3_OBJS)
	$(CXX) -o $@ $(TEST3_OBJS)

$(TEST4) : $(TEST4_OBJS)
	$(CXX) -o $@ $(TEST4_OBJS)

$(TEST5) : $(TEST5_OBJS)
	$(CXX) -o $@ $(TEST5_OBJS)

new :
	$(MAKE) clean
	$(MAKE)
# DO NOT DELETE

./src/main.o: src/include/mcm.hpp src/include/fourier_path.hpp
./src/main.o: src/include/fourier.hpp src/include/fourier/fourier.tpp
./src/main.o: src/include/fourier_path/fourier_path.tpp
./src/main.o: src/include/action.hpp src/include/node_weight_table.hpp
./src/main.o: src/include/action/action.tpp src/include/mcm/mcm_basic.tpp
./src/main.o: src/include/mcm/mcm_move.tpp src/include/mcm/mcm_opt.tpp
./src/test/test_action_kepler.o: src/include/action.hpp
./src/test/test_action_kepler.o: src/include/fourier_path.hpp
./src/test/test_action_kepler.o: src/include/fourier.hpp
./src/test/test_action_kepler.o: src/include/fourier/fourier.tpp
./src/test/test_action_kepler.o: src/include/fourier_path/fourier_path.tpp
./src/test/test_action_kepler.o: src/include/node_weight_table.hpp
./src/test/test_action_kepler.o: src/include/action/action.tpp
./src/test/test_action_simple.o: src/include/action.hpp
./src/test/test_action_simple.o: src/include/fourier_path.hpp
./src/test/test_action_simple.o: src/include/fourier.hpp
./src/test/test_action_simple.o: src/include/fourier/fourier.tpp
./src/test/test_action_simple.o: src/include/fourier_path/fourier_path.tpp
./src/test/test_action_simple.o: src/include/node_weight_table.hpp
./src/test/test_action_simple.o: src/include/action/action.tpp
./src/test/test_action_vaildity.o: src/include/action.hpp
./src/test/test_action_vaildity.o: src/include/fourier_path.hpp
./src/test/test_action_vaildity.o: src/include/fourier.hpp
./src/test/test_action_vaildity.o: src/include/fourier/fourier.tpp
./src/test/test_action_vaildity.o: src/include/fourier_path/fourier_path.tpp
./src/test/test_action_vaildity.o: src/include/node_weight_table.hpp
./src/test/test_action_vaildity.o: src/include/action/action.tpp
./src/test/test_fourier.o: src/include/fourier.hpp
./src/test/test_fourier.o: src/include/fourier/fourier.tpp
./src/test/test_fourier_path.o: src/include/fourier_path.hpp
./src/test/test_fourier_path.o: src/include/fourier.hpp
./src/test/test_fourier_path.o: src/include/fourier/fourier.tpp
./src/test/test_fourier_path.o: src/include/fourier_path/fourier_path.tpp
