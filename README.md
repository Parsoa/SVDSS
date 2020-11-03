# Sample-specific string detection from accurate long reads
Efficient computation of A-specific string w.r.t. a set {B,C,...,Z} of other long reads samples. A A-specific string is a string which occur only in sample A and not in the others. 

##### Use-Cases
* compute strings specific to child w.r.t. parents
* compute strings specific to individual A from population P<sub>A</sub> w.r.t. individual B from population P<sub>B</sub>

## Dependencies

C++11-compliant compiler and [ropebwt2](https://github.com/lh3/ropebwt2) library. For convenience, ropebwt2 is included in the repository.

## Download and Installation
```
git clone --recursive https://github.com/Parsoa/Variation-Discovery.git
cd Variation-Discovery/ropebwt2
make
cd ..
make
```

## How-To
Let's assume to have 3 samples: A, B, and, C. If we want to compute A-specific strings we have to:

1. index samples B and C:
```
./stella pingpong index -b /path/to/sample/B > B.index.bin
./stella pingpong -a B.index.bin /path/to/sample/C > BC.index.fmd
```
2. search A-specific strings
```
./stella pingpong search [B.index.bin] /path/to/sample/A [nthreads]
```
3. the search step splits the input sample in batches, to combine them in a single fastq file:
```
./aggregate /path/to/curr/dir > specific.fastq
```

### Usage
```
Usage: stella pingpong index [-h] [-b] [-a index] <sample>

Optional arguments:
      -h, --help            display this help and exit
      -b, --binary          output index in binary format
      -a, --append          append to existing index (must be stored in binary)
      -t, --threads         number of threads (default:1)

Usage: stella pingpong search [-h] <index> <sample> <threads>
```

##### Notes
* to append (`-a`) to an existing index, the existing index must be stored in binary format (`-b` option)
* to query the index, it must be stored in FMD format (default)
* the search output is stored in multiple `solution_batch_*.fastq` files (created in the current directory)

### Example
```
./stella pingpong index -b example/father.fq > example/father.fq.bin
./stella pingpong index -a example/father.fq.bin example/mother.fq > example/index.fmd
./stella pingpong search example/index.fmd example/child.fq 1
 ```
 
 ### Authors
 For inquiries on this software please open an [issue](https://github.com/Parsoa/Variation-Discovery/issues) or contact either [Parsoa Khorsand](https://github.com/parsoa) or [Luca Denti](https://github.com/ldenti/).
