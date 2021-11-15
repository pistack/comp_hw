include make.inc

.SUFFIXES: .cpp .o

INC = -I include/

DFLAGS = -D MONITOR=$(MONITOR) -D PRECISION_LEVEL=$(PRECISION_LEVEL) -D $(PATH_TYPE)
CPPFLAGS = $(DFLAGS) $(INC)

OBJS = ./src/main.o 
SRCS = ./src/main.cpp

TARGET = ./bin/hw3

export DFLAGS

all : $(TARGET)

$(TARGET) : $(OBJS)
	$(CXX) -o $@ $(OBJS)

dep :
	makedepend $(INC) $(SRCS)
	$(MAKE) -C test dep

clean :
	rm -rf $(OBJS) $(TARGET)
	$(MAKE) -C test clean

.PHONY : test
test : 
	$(MAKE) -C test test

new :
	$(MAKE) clean
	$(MAKE)
# DO NOT DELETE

./src/main.o: include/libpath/math_const.hpp
./src/main.o: include/libpath/math_const/math_const.tpp
./src/main.o: include/libmcm/mcm_fourier.hpp include/libpath/fourier_path.hpp
./src/main.o: include/libpath/fourier.hpp include/libpath/math_const.hpp
./src/main.o: include/libpath/fourier/fourier.tpp
./src/main.o: include/libpath/fourier_path/fourier_path.tpp
./src/main.o: include/libmcm/mcm.hpp include/libpath/action.hpp
./src/main.o: include/libpath/node_weight_table.hpp
./src/main.o: include/libpath/action/action.tpp
./src/main.o: include/libmcm/mcm/mcm_basic.tpp
./src/main.o: include/libmcm/mcm/mcm_move.tpp include/libmcm/mcm/mcm_opt.tpp
./src/main.o: include/libmcm/mcm_fourier/mcm_fourier.tpp
