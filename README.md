# Sample-specific string detection from accurate long reads
Efficient computation of A-specific string w.r.t. a set {B,C,...,Z} of other long reads samples. A A-specific string is a string which occur only in sample A and not in the others. 

##### Use-Cases
* compute strings specific to child w.r.t. parents
* compute strings specific to individual A from population P<sub>A</sub> w.r.t. individual B from population P<sub>B</sub>

## Dependencies

C++11-compliant compiler and [ropebwt2](https://github.com/lh3/ropebwt2) library. For convenience, ropebwt2 is included in the repository.

## Download and Installation
```
git clone --recursive https://github.com/Parsoa/PingPong.git
cd PingPong/ropebwt2
make
cd ..
make
```

## How-To
Let's assume we have 3 samples A, B, and, C. To compute A-specific strings we have to:

1. index samples B and C:
```
./stella pingpong index -b /path/to/sample/B > B.index.bin
./stella pingpong index -a B.index.bin /path/to/sample/C > BC.index.fmd
```
2. search for A-specific strings in the index
```
./stella pingpong search --index [B.index.bin] --fastq /path/to/sample/A --threads [nthreads]
```
3. the search step splits the input sample in batches, to aggregate them to a single fastq file:
```
./stella aggregate --workdir /path/to/string/batches --threads
```

The algorithm will out a file named `solution.aggregated.fastq` with the list of A-specific strings that are repeated more than

### PingPong Algorithm Usage
```
Usage: stella pingpong index [--binary] [--append index] --fastq /path/to/fastq/file [--threads threads]

Optional arguments:
      -b, --binary          output index in binary format
      -a, --append          append to existing index (must be stored in binary)
      -t, --threads         number of threads (default:1)

Usage: stella pingpong search [--index index] [--fastq fastq] [--threads threads]
```

##### Notes
* to append (`-a`) to an existing index, the existing index must be stored in binary format (`-b` option)
* to query the index, it must be stored in FMD format (default)
* the search output is stored in multiple `solution_batch_*.fastq` files (created in the current directory)

### Example

```
./stella pingpong index --binary --fastq example/father.fq > example/father.fq.bin
./stella pingpong index --append --fastq example/father.fq.bin example/mother.fq > example/index.fmd
./stella pingpong search --index example/index.fmd --fastq example/child.fq --threads 4
 ```
 
 ### Authors
 For inquiries on this software please open an [issue](https://github.com/Parsoa/Variation-Discovery/issues) or contact either [Parsoa Khorsand](https://github.com/parsoa) or [Luca Denti](https://github.com/ldenti/).
