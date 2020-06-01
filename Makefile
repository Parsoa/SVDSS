ROPE_PREFIX=./ropebwt2
CXX=g++
CXXFLAGS=-std=c++11 -Wall -O3 -I$(ROPE_PREFIX) -fopenmp -Wno-comment
LIBS=-lm -lz -lpthread
OBJS=$(ROPE_PREFIX)/mrope.o $(ROPE_PREFIX)/rope.o $(ROPE_PREFIX)/rld0.o $(ROPE_PREFIX)/rle.o

all: main

main: main.cpp
	$(CXX) -g $(CXXFLAGS) $(OBJS) $< -o main $(LIBS) 

debug: main.cpp
	$(CXX) -g $(CXXFLAGS) $(OBJS) -DDEBUG_MODE $< -o debug $(LIBS) 

clean:
	rm -rf main debug
