CXX=g++
CXXFLAGS=-std=c++14 -Wall -O3 -I/pasteur/sonic/homes/ldenti/software/abPOA/include -I/pasteur/sonic/homes/ldenti/software/parasail -I/pasteur/sonic/homes/ldenti/code/PingPong/htslib -I/pasteur/sonic/homes/ldenti/software/interval-tree
LIBS=-L/pasteur/sonic/homes/ldenti/software/abPOA/lib -L/pasteur/sonic/homes/ldenti/software/parasail/build -L/pasteur/sonic/homes/ldenti/code/PingPong/htslib
LDFLAGS=-lz -labpoa -lparasail -lhts

all: main

main: sv.o clipler.o main.o
	@echo "* Linking $@"
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o main
