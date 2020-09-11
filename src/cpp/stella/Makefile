CXX = g++
LDFLAGS = -lm -lpthread
CXXFLAGS = -std=c++17 -O3 -funroll-loops -fopenmp

HPP = $(wildcard *.hpp)
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

all: stella

debug: CXXFLAGS += -DDEBUG -g
debug: stella 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

stella: $(OBJ) $(HPP) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDFLAGS)

clean:
	rm *.o
.PHONY: clean

