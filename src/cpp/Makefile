CC = g++
CFLAGS = -O3 -funroll-loops -std=c++17 -I htslib -fconcepts
LDFLAGS = /software/htslib/1.8/lssc0-linux/lib/libhts.a -lz -lm -lpthread -llzma -lbz2 -lcurl -lcrypto -DSEQAN3_HAS_ZLIB=1 -DSEQAN3_HAS_BZIP2=1 -lstdc++fs -lsdsl -ldivsufsort -ldivsufsort64 -L /home/pkhorsand/lib
INCFLAGS = -I /home/pkhorsand/include
# Docker config
DCFLAGS = -O3 -funroll-loops -std=c++11
DLDFLAGS = -lz -lm -lpthread -lhts
SOURCES = counter.cpp
OBJECTS = $(SOURCES:.cpp=.o)

counter.o: counter.cpp
	$(CC) -g $(CFLAGS) $(INCFLAGS) -c $^ -o $@

counter: counter.o
	$(CC) counter.o -o counter.out $(LDFLAGS)

search.o: search.cpp
	$(CC) -g $(CFLAGS) $(INCFLAGS) -c $^ -o $@

search: search.o
	$(CC) search.o -o search.out $(LDFLAGS)

debug: counter.cpp
	$(CC) -g $(CFLAGS) -DDEBUG counter.cpp -o counter.out $(LDFLAGS)

docker: counter.cpp
	$(CC) -g $(DCFLAGS) counter.cpp -o counter.out $(DLDFLAGS)

