CXX=g++
CXXFLAGS=-std=c++11 -Wall -Wno-unused-variable -O3 -I../abPOA/include -I../lib/parasail/include -I../htslib -I../lib/intervaltree
LIBS=-L../abPOA/lib -L../lib/parasail/lib -L../htslib
LDFLAGS=-lz -labpoa -lparasail -lhts

all: main

main: main.o
	@echo "* Linking $@"
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o main
