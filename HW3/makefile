.SUFFIXES: .cpp .o

CXX = g++
INC = -I include/
# CPPFLAGS = -O -g -Wall -Werror -pedantic -D MONITOR=0 -D PRECISION_LEVEL=1 $(INC)
CPPFLAGS = -O3 -D MONITOR=1 -D PRECISION_LEVEL=0 $(INC) # for single precision
# CPPFLAGS = -O3 -D MONITOR=1 -D PRECISION_LEVEL=1 $(INC) # for double precision
# CPPFLAGS = -O3 -D MONITOR=0 -D PRECISION_LEVEL=0 $(INC) # for single precision
# CPPFLAGS = -O3 -D MONITOR=0 -D PRECISION_LEVEL=1 $(INC) # for double precision

OBJS = ./src/main.o 
TEST1_OBJS = ./test/test_action_vaildity.o
TEST2_OBJS = ./test/test_action_simple.o
TEST3_OBJS = ./test/test_action_kepler.o 
TEST4_OBJS = ./test/test_fourier.o 
TEST5_OBJS = ./test/test_fourier_path.o
TEST6_OBJS = ./test/test_bezier.o
TEST7_OBJS = ./test/test_bezier_path.o

TEST_OBJS = $(TEST1_OBJS) $(TEST2_OBJS) $(TEST3_OBJS) $(TEST4_OBJS) \
$(TEST5_OBJS) $(TEST6_OBJS) $(TEST7_OBJS)

SRCS = ./src/main.cpp ./test/test_action_kepler.cpp \
./test/test_action_simple.cpp ./test/test_action_vaildity.cpp \
./test/test_fourier.cpp ./test/test_fourier_path.cpp \
./test/test_bezier.cpp ./test/test_bezier_path.cpp

TARGET = ./bin/hw3
TEST1 = test_action_vaildity
TEST2 = test_action_simple
TEST3 = test_action_kepler
TEST4 = test_fourier
TEST5 = test_fourier_path
TEST6 = test_bezier
TEST7 = test_bezier_path 

TEST = $(TEST1) $(TEST2) $(TEST3) $(TEST4) $(TEST5) \
$(TEST6) $(TEST7)


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
	./$(TEST6)
	./$(TEST7)
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

$(TEST6) : $(TEST6_OBJS)
	$(CXX) -o $@ $(TEST6_OBJS)

$(TEST7) : $(TEST7_OBJS)
	$(CXX) -o $@ $(TEST7_OBJS)

new :
	$(MAKE) clean
	$(MAKE)
# DO NOT DELETE

./src/main.o: include/libmcm/mcm.hpp include/libpath/fourier_path.hpp
./src/main.o: include/libpath/fourier.hpp include/libpath/fourier/fourier.tpp
./src/main.o: include/libpath/fourier_path/fourier_path.tpp
./src/main.o: include/libpath/action.hpp
./src/main.o: include/libpath/node_weight_table.hpp
./src/main.o: include/libpath/action/action.tpp
./src/main.o: include/libmcm/mcm/mcm_basic.tpp
./src/main.o: include/libmcm/mcm/mcm_move.tpp include/libmcm/mcm/mcm_opt.tpp
./test/test_action_kepler.o: include/libpath/fourier_path.hpp
./test/test_action_kepler.o: include/libpath/fourier.hpp
./test/test_action_kepler.o: include/libpath/fourier/fourier.tpp
./test/test_action_kepler.o: include/libpath/fourier_path/fourier_path.tpp
./test/test_action_kepler.o: include/libpath/action.hpp
./test/test_action_kepler.o: include/libpath/node_weight_table.hpp
./test/test_action_kepler.o: include/libpath/action/action.tpp
./test/test_action_simple.o: include/libpath/bezier_path.hpp
./test/test_action_simple.o: include/libpath/bezier.hpp
./test/test_action_simple.o: include/libpath/bezier/bezier.tpp
./test/test_action_simple.o: include/libpath/bezier_path/bezier_path.tpp
./test/test_action_simple.o: include/libpath/fourier_path.hpp
./test/test_action_simple.o: include/libpath/fourier.hpp
./test/test_action_simple.o: include/libpath/fourier/fourier.tpp
./test/test_action_simple.o: include/libpath/fourier_path/fourier_path.tpp
./test/test_action_simple.o: include/libpath/action.hpp
./test/test_action_simple.o: include/libpath/node_weight_table.hpp
./test/test_action_simple.o: include/libpath/action/action.tpp
./test/test_action_vaildity.o: include/libpath/fourier_path.hpp
./test/test_action_vaildity.o: include/libpath/fourier.hpp
./test/test_action_vaildity.o: include/libpath/fourier/fourier.tpp
./test/test_action_vaildity.o: include/libpath/fourier_path/fourier_path.tpp
./test/test_action_vaildity.o: include/libpath/action.hpp
./test/test_action_vaildity.o: include/libpath/node_weight_table.hpp
./test/test_action_vaildity.o: include/libpath/action/action.tpp
./test/test_fourier.o: include/libpath/fourier.hpp
./test/test_fourier.o: include/libpath/fourier/fourier.tpp
./test/test_fourier_path.o: include/libpath/fourier_path.hpp
./test/test_fourier_path.o: include/libpath/fourier.hpp
./test/test_fourier_path.o: include/libpath/fourier/fourier.tpp
./test/test_fourier_path.o: include/libpath/fourier_path/fourier_path.tpp
./test/test_bezier.o: include/libpath/bezier.hpp
./test/test_bezier.o: include/libpath/bezier/bezier.tpp
./test/test_bezier_path.o: include/libpath/bezier_path.hpp
./test/test_bezier_path.o: include/libpath/bezier.hpp
./test/test_bezier_path.o: include/libpath/bezier/bezier.tpp
./test/test_bezier_path.o: include/libpath/bezier_path/bezier_path.tpp
