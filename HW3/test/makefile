.SUFFIXES: .cpp .o


INC = -I ../include/ -I include/

# DEBUG
CPPFLAGS = -O -g -Wall -pedantic $(DFLAGS) $(INC)

TEST1_OBJS = test_action_vaildity.o
TEST2_OBJS = test_action_simple.o
TEST3_OBJS = test_action_kepler.o 
TEST4_OBJS = test_fourier.o 
TEST5_OBJS = test_fourier_path.o
TEST6_OBJS = test_bezier.o
TEST7_OBJS = test_bezier_path.o

TEST_OBJS = test.o $(TEST1_OBJS) $(TEST2_OBJS) $(TEST3_OBJS) $(TEST4_OBJS) \
$(TEST5_OBJS) $(TEST6_OBJS) $(TEST7_OBJS)

SRCS = test.cpp test_action_kepler.cpp \
test_action_simple.cpp test_action_vaildity.cpp \
test_fourier.cpp test_fourier_path.cpp \
test_bezier.cpp test_bezier_path.cpp

TARGET = build test


all : $(TARGET)

build : $(TEST_OBJS)
	$(CXX) -o test $(TEST_OBJS)

test: build
	./test

dep :
	makedepend $(INC) $(SRCS)

clean :
	rm -rf test $(TEST_OBJS)

new :
	$(MAKE) clean
	$(MAKE)

# DO NOT DELETE

test.o: include/test.hpp
test_action_kepler.o: ../include/libpath/fourier_path.hpp
test_action_kepler.o: ../include/libpath/fourier.hpp
test_action_kepler.o: ../include/libpath/fourier/fourier.tpp
test_action_kepler.o: ../include/libpath/fourier_path/fourier_path.tpp
test_action_kepler.o: ../include/libpath/action.hpp
test_action_kepler.o: ../include/libpath/node_weight_table.hpp
test_action_kepler.o: ../include/libpath/action/action.tpp include/test.hpp
test_action_simple.o: ../include/libpath/bezier_path.hpp
test_action_simple.o: ../include/libpath/bezier.hpp
test_action_simple.o: ../include/libpath/bezier/bezier.tpp
test_action_simple.o: ../include/libpath/bezier_path/bezier_path.tpp
test_action_simple.o: ../include/libpath/fourier_path.hpp
test_action_simple.o: ../include/libpath/fourier.hpp
test_action_simple.o: ../include/libpath/fourier/fourier.tpp
test_action_simple.o: ../include/libpath/fourier_path/fourier_path.tpp
test_action_simple.o: ../include/libpath/action.hpp
test_action_simple.o: ../include/libpath/node_weight_table.hpp
test_action_simple.o: ../include/libpath/action/action.tpp include/test.hpp
test_action_vaildity.o: ../include/libpath/fourier_path.hpp
test_action_vaildity.o: ../include/libpath/fourier.hpp
test_action_vaildity.o: ../include/libpath/fourier/fourier.tpp
test_action_vaildity.o: ../include/libpath/fourier_path/fourier_path.tpp
test_action_vaildity.o: ../include/libpath/action.hpp
test_action_vaildity.o: ../include/libpath/node_weight_table.hpp
test_action_vaildity.o: ../include/libpath/action/action.tpp include/test.hpp
test_fourier.o: ../include/libpath/fourier.hpp
test_fourier.o: ../include/libpath/fourier/fourier.tpp include/test.hpp
test_fourier_path.o: ../include/libpath/fourier_path.hpp
test_fourier_path.o: ../include/libpath/fourier.hpp
test_fourier_path.o: ../include/libpath/fourier/fourier.tpp
test_fourier_path.o: ../include/libpath/fourier_path/fourier_path.tpp
test_fourier_path.o: include/test.hpp
test_bezier.o: ../include/libpath/bezier.hpp
test_bezier.o: ../include/libpath/bezier/bezier.tpp include/test.hpp
test_bezier_path.o: ../include/libpath/bezier_path.hpp
test_bezier_path.o: ../include/libpath/bezier.hpp
test_bezier_path.o: ../include/libpath/bezier/bezier.tpp
test_bezier_path.o: ../include/libpath/bezier_path/bezier_path.tpp
test_bezier_path.o: include/test.hpp
