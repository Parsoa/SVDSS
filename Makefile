ROPE_PREFIX=./ropebwt2
CXX=g++
CXXFLAGS=-std=c++17 -Wall -O3 -I$(ROPE_PREFIX) -fopenmp -Wno-comment
LDFLAGS= -lm -lpthread -lz
OBJS=$(ROPE_PREFIX)/mrope.o $(ROPE_PREFIX)/rope.o $(ROPE_PREFIX)/rld0.o $(ROPE_PREFIX)/rle.o

HPP = $(wildcard *.hpp)
SRC = $(wildcard *.cpp)
OBJS := $(OBJS) $(SRC:.cpp=.o)

all: stella 
debug: CXXFLAGS += -DDEBUG -g
debug: stella 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

stella: $(OBJS) $(HPP) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm *.o
.PHONY: clean

