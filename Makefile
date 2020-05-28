CXX=g++
CXXFLAGS=-std=c++11 -Wall -O3 -I./fermi -fopenmp -Wno-comment
LIBS=-lm -lz -lpthread
OBJS=fermi/utils.o fermi/seq.o fermi/ksa.o fermi/ksa64.o fermi/rld.o fermi/exact.o fermi/merge.o fermi/sub.o fermi/correct.o \
	fermi/build.o fermi/smem.o fermi/unitig.o fermi/seqsort.o fermi/cmp.o fermi/cmd.o fermi/example.o \
	fermi/ksw.o fermi/mag.o fermi/bubble.o fermi/scaf.o fermi/bcr.o fermi/bprope6.o fermi/ropebwt.o \

all: main

luca: luca.cpp
	$(CXX) -g $(CXXFLAGS) -I./KMC $(OBJS) -DNDEBUG KMC/kmc_api/kmc_file.o KMC/kmc_api/kmer_api.o KMC/kmc_api/mmer.o $< -o luca.out $(LIBS) 

main: main.cpp
	$(CXX) -g $(CXXFLAGS) $(OBJS) $< -o main.out $(LIBS) 

debug: main.cpp
	$(CXX) -g $(CXXFLAGS) $(OBJS) -DDEBUG_MODE $< -o debug.out $(LIBS) 

#main: main.o
#	@echo "* Linking $@"
#	$(CXX) $(CXXFLAGS) $(OBJS) $^ -o $@ $(LIBS)
#
#%.o: %.cpp
#	@echo '* Compiling $<'
#	$(CXX) $(CXXFLAGS) -o $@ -c $<
#
#clean:
#	rm -rf main.o
