KMC_PREFIX=./KMC

CXX=g++
CXXFLAGS=-std=c++14 -Wall -O3 -I. -I$(KMC_PREFIX) -Wno-char-subscripts
LIBS=-lz
OBJS=$(KMC_PREFIX)/kmc_api/kmc_file.o $(KMC_PREFIX)/kmc_api/kmer_api.o $(KMC_PREFIX)/kmc_api/mmer.o

.PHONY: all

all: cutsample

cutsample: main.o
	@echo "* Linking $@"
	$(CXX) $(CXXFLAGS) $(OBJS) $^ -o $@ $(LIBS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf cutsample main.o
