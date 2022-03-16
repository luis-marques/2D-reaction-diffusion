CC			= g++
CXXFLAGS	= -std=c++11 -Wall -Wshadow -O0 -g
LDLIBS		= -llapack -lblas -lboost_program_options
TARGET		= main
HDRS 		= ReactionDiffusion.h

TEST1_PARAMS = --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.75 --b 0.06 --eps 50.0 --mu1 5.0 --mu2 0.0
TEST2_PARAMS = --dt 0.001 --T 100 --Nx 251 --Ny 251 --a 0.75 --b 0.06 --eps 13.0 --mu1 5.0 --mu2 0.0
TEST3_PARAMS = --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.5 --b 0.1 --eps 50.0 --mu1 5.0 --mu2 0.0
TEST4_PARAMS = --dt 0.001 --T 100 --Nx 151 --Ny 81 --a 0.75 --b 0.0001 --eps 12.5 --mu1 1.0 --mu2 0.01


default: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).o ReactionDiffusion.o

%.o: %.cpp $(HDRS)

.PHONY: clean run

run: $(TARGET)
	./$(TARGET) --dt 0.001 --T 100

clean:
	rm -f $(TARGET) *.o
