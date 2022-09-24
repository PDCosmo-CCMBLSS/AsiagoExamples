# Set compiler
CC = gcc
CPPFLAGS = -lstdc++ -std=c++17 -g -fPIC

#=======================================================
# Includes and libraries
#=======================================================
INC = -I/usr/local/include -Iinclude/ -I./
LIBS = -L/usr/local/lib -lm

#-------------------------------------------------------
# Boost
INC += -I/usr/local/lib/boost_1_68_0/ -I../../../lib/boost_1_68_0/

#=======================================================
# Options
#=======================================================
OPTIONS = -Wall

# Optimization
OPTIONS = # -O3 -march=native
OPTIONS += -O3

#=======================================================
# Object files
#=======================================================
OBJS = background.o




comdist: main.cpp $(OBJS)
	${CC} -o $@ $^ $(CPPFLAGS) $(OPTIONS) $(INC) $(LIBS)

comdist_best_practice: main_best_practice.cpp $(OBJS)
	${CC} -o $@ $^ $(CPPFLAGS) $(OPTIONS) $(INC) $(LIBS)
	
comdist_error: main_error.cpp $(OBJS)
	${CC} -o $@ $^ $(CPPFLAGS) $(OPTIONS) $(INC) $(LIBS)

%.o: %.cpp
	${CC} -c -o $@ $< $(CPPFLAGS) $(OPTIONS) $(INC)
	
clean:
	rm -rf ./*.o comdist comdist_best_practice
