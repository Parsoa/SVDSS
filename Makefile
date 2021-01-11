ROPE_PREFIX=./ropebwt2
HTS_PREFIX=./htslib
CXX=g++
CXXFLAGS=-std=c++11 -Wall -Wno-comment -O3 -I$(ROPE_PREFIX) -I$(HTS_PREFIX) -fopenmp -Wno-sign-compare -Wno-unused-variable
LIBS=-L$(HTS_PREFIX) -lm -lpthread -lz -lhts
BIN=$(ROPE_PREFIX)/mrope.o $(ROPE_PREFIX)/rope.o $(ROPE_PREFIX)/rld0.o $(ROPE_PREFIX)/rle.o

HPP = $(wildcard *.hpp)
SRC = $(wildcard *.cpp)
OBJS = $(BIN) $(SRC:.cpp=.o)

LOCAL_BUILD ?= 0
ifneq (${LOCAL_BUILD}), 0)
	HPP += $(wildcard src/*.hpp)
	SRC += $(wildcard src/*.cpp)
	CXX += -Isrc -DLOCAL_BUILD
endif	

$(info $(LOCAL_BUILD))

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

