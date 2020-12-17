ROPE_PREFIX=./ropebwt2
HTSLIBFLAGS =-lhts -lz -llzma -lbz2 -lcurl -lcrypto
CXX=g++-10
CXXFLAGS=-std=c++11 -Wall -O3 -I$(ROPE_PREFIX) -I htslib -fopenmp -Wno-comment
LDFLAGS= -lm -lpthread -lz -L/usr/local/lib -L/usr/local/opt/openssl/lib $(HTSLIBFLAGS)
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

