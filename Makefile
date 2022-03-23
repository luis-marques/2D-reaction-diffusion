# The "DEBUG" variable is to facilitate changing between "debug"
# and "release" builds of the code. The first containts debugging symbols,
# compiler flags to provide warnings and is compiled without any optimizations.
# The second is the default build of the code and contains optimizations flags.
# To compile the "debug" build run: "$ make DEBUG=1".
DEBUG			?= 0
ifeq ($(DEBUG), 1)
	CXXFLAGS	= -fopenmp -std=c++11 -Wall -Wshadow -pedantic -O0 -g
else
	CXXFLAGS	= -fopenmp -std=c++11 -O3 -march=native
endif

CC				= g++
LDLIBS			= -lboost_program_options
LDFLAGS			= -fopenmp
TARGET			= main
HDRS 			= ReactionDiffusion.h
DOXY 			= doxygen
DOXYFLAGS		= Doxyfile

# Command-line arguments providing parameters of the 4 test cases
TEST1_PARAMS 	= --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.75 --b 0.06 --eps 50.0 --mu1 5.0 --mu2 0.0
TEST2_PARAMS 	= --dt 0.001 --T 100 --Nx 251 --Ny 251 --a 0.75 --b 0.06 --eps 13.0 --mu1 5.0 --mu2 0.0
TEST3_PARAMS 	= --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.5 --b 0.1 --eps 50.0 --mu1 5.0 --mu2 0.0
TEST4_PARAMS 	= --dt 0.001 --T 100 --Nx 151 --Ny 81 --a 0.75 --b 0.0001 --eps 12.5 --mu1 1.0 --mu2 0.01

# Setting the "default" target
default: all
all: $(TARGET) doc

$(TARGET): $(TARGET).o ReactionDiffusion.o

# Header file change requires re-compilation (both .cpps include the header file)
%.o: %.cpp $(HDRS)

# These make targets don't create new files but rather run commands
.PHONY: clean run doc test1 test2 test3 test4 debug

run: $(TARGET)
	./$(TARGET) --dt 0.001 --T 100

clean:
	rm -f $(TARGET) *.o
	rm -rf html
	rm -rf latex

doc: *.cpp
	$(DOXY) $(DOXYFILE)

test1: $(TARGET)
	./$(TARGET) $(TEST1_PARAMS)

test2: $(TARGET)
	./$(TARGET) $(TEST2_PARAMS)
	
test3: $(TARGET)
	./$(TARGET) $(TEST3_PARAMS)
	
test4: $(TARGET)
	./$(TARGET) $(TEST4_PARAMS)