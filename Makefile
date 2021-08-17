ROPE_PREFIX=./ropebwt2
HTS_PREFIX=./htslib
RF_PREFIX=./rapidfuzz-cpp
ABPOA_PREFIX=./abPOA
KSW_PREFIX=./ksw2
INTTREE_PREFIX=./interval-tree
CXX=g++
CXXFLAGS=-std=c++14 -Wall -Wno-comment -O3 -I$(ROPE_PREFIX) -I$(HTS_PREFIX) -I$(RF_PREFIX) -I$(ABPOA_PREFIX)/include -I$(KSW_PREFIX) -I$(INTTREE_PREFIX) -fopenmp -Wno-sign-compare -Wno-unused-variable
LIBS=-L$(HTS_PREFIX) -L$(ABPOA_PREFIX)/lib -lm -lpthread -lz -lhts -labpoa
BIN=$(ROPE_PREFIX)/mrope.o $(ROPE_PREFIX)/rope.o $(ROPE_PREFIX)/rld0.o $(ROPE_PREFIX)/rle.o
KSW_OBJS=$(addprefix $(KSW_PREFIX)/, kalloc.o ksw2_gg.o ksw2_gg2.o ksw2_gg2_sse.o ksw2_extz.o ksw2_extz2_sse.o ksw2_extd.o ksw2_extd2_sse.o ksw2_extf2_sse.o ksw2_exts2_sse.o)

HPP = $(wildcard *.hpp)
SRC = $(wildcard *.cpp)
OBJS = $(BIN) $(KSW_OBJS) $(SRC:.cpp=.o)

LOCAL_BUILD ?= 0
ifneq ($(LOCAL_BUILD), 0)
	HPP += $(wildcard src/*.hpp)
	SRC += $(wildcard src/*.cpp)
	CXX += -Isrc -DLOCAL_BUILD
endif	

$(info $(LOCAL_BUILD))

all: PingPong 
debug: CXXFLAGS += -DDEBUG_MODE -g
debug: PingPong 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

PingPong: $(OBJS) $(HPP) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f PingPong *.o
.PHONY: clean

