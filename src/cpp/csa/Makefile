CC = g++
CFLAGS = -O3 -funroll-loops -std=c++17 -fconcepts -DSEQAN3_HAS_ZLIB=1 -DSEQAN3_HAS_BZIP2=1 
LDFLAGS = /software/htslib/1.8/lssc0-linux/lib/libhts.a -lz -lm -lpthread -llzma -lbz2 -lcurl -lcrypto -lstdc++fs -lsdsl -ldivsufsort -ldivsufsort64 -L /home/pkhorsand/lib
INCFLAGS = -I /home/pkhorsand/include -I htslib 
# Docker config
SOURCES = counter.cpp
OBJECTS = $(SOURCES:.cpp=.o)

counter.o: counter.cpp
	$(CC) -g $(CFLAGS) $(INCFLAGS) -c $^ -o $@

counter: counter.o
	$(CC) counter.o -o counter.out $(LDFLAGS)

debug: counter.cpp
	$(CC) -g $(CFLAGS) $(INCFLAGS) -DDEBUG_BUILD counter.cpp -o counter.out $(LDFLAGS)
