ROPE_PREFIX=./ropebwt2
HTSLIBFLAGS =/software/htslib/1.8/lssc0-linux/lib/libhts.a -lz -llzma -lbz2 -lcurl -lcrypto
CXX=g++
CXXFLAGS=-std=c++11 -Wall -O3 -I$(ROPE_PREFIX) -I htslib -fopenmp -Wno-comment
LDFLAGS= -lm -lpthread -lz $(HTSLIBFLAGS)
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
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm *.o
.PHONY: clean

