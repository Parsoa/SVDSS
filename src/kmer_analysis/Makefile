CXXFLAGS=-DNDEBUG -Wall -O3 -std=c++11 -I. -I./include -I./KMC -fopenmp
LIBS=-L./lib -lz -lsdsl -ldivsufsort -ldivsufsort64

.PHONY: all

all: main

main: main.o ./KMC/kmc_api/kmc_file.o ./KMC/kmc_api/kmer_api.o ./KMC/kmc_api/mmer.o
	@echo "* Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o
