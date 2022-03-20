# For ease of changing between "debug" and "release" builds
# Default behavior is optimized build, to get debug build run "$ make DEBUG=1"
DEBUG		?= 0

# https://wiki.gentoo.org/wiki/GCC_optimization#Optimizing

ifeq ($(DEBUG), 1)
	CXXFLAGS	= -fopenmp -std=c++11 -Wall -Wshadow -pedantic -O0 -g
else
	CXXFLAGS	= -fopenmp -std=c++11 -O4 -march=native
endif

CC			= g++
LDLIBS		= -lboost_program_options
LDFLAGS		= -fopenmp
TARGET		= main
HDRS 		= ReactionDiffusion.h

# Command-line arguments providing parameters of the 4 test cases
TEST1_PARAMS = --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.75 --b 0.06 --eps 50.0 --mu1 5.0 --mu2 0.0
TEST2_PARAMS = --dt 0.001 --T 100 --Nx 251 --Ny 251 --a 0.75 --b 0.06 --eps 13.0 --mu1 5.0 --mu2 0.0
TEST3_PARAMS = --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.5 --b 0.1 --eps 50.0 --mu1 5.0 --mu2 0.0
TEST4_PARAMS = --dt 0.001 --T 100 --Nx 151 --Ny 81 --a 0.75 --b 0.0001 --eps 12.5 --mu1 1.0 --mu2 0.01

# Setting "default" target
default: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).o ReactionDiffusion.o

# Header file change requires re-compilation
%.o: %.cpp $(HDRS)

# These make targets don't create new files but rather run commands
.PHONY: clean run test1 test2 test3 test4 debug

run: $(TARGET)
	./$(TARGET) --dt 0.001 --T 100

clean:
	rm -f $(TARGET) *.o

test1: $(TARGET)
	./$(TARGET) $(TEST1_PARAMS)

test2: $(TARGET)
	./$(TARGET) $(TEST2_PARAMS)
	
test3: $(TARGET)
	./$(TARGET) $(TEST3_PARAMS)
	
test4: $(TARGET)
	./$(TARGET) $(TEST4_PARAMS)