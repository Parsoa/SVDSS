ROPE_PREFIX=./ropebwt2
HTS_PREFIX=./htslib
CXX=g++
CXXFLAGS=-std=c++11 -Wall -Wno-comment -O3 -I$(ROPE_PREFIX) -I$(HTS_PREFIX) -fopenmp
LIBS=-L$(HTS_PREFIX) -lm -lpthread -lz -lhts
OBJS=$(ROPE_PREFIX)/mrope.o $(ROPE_PREFIX)/rope.o $(ROPE_PREFIX)/rld0.o $(ROPE_PREFIX)/rle.o

HPP = $(wildcard *.hpp)
SRC = $(wildcard *.cpp)
OBJS := $(OBJS) $(SRC:.cpp=.o)

all: stella 
debug: CXXFLAGS += -DDEBUG_MODE -g
debug: stella 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

stella: $(OBJS) $(HPP) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f stella *.o
.PHONY: clean

