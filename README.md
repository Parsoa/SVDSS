![C/C++ CI](https://github.com/Parsoa/PingPong/workflows/C/C++%20CI/badge.svg)

# Sample-specific string detection from accurate long reads
Efficient computation of A-specific string w.r.t. a set {B,C,...,Z} of other long reads samples. A A-specific string is a string which occur only in sample A and not in the others. 

**Note:** if you are looking for the _biorxiv_ implementation, please refer to [this release](https://github.com/Parsoa/PingPong/releases/tag/v1.0.0-pingpong).

##### Use-Cases
* compute strings specific to child w.r.t. parents
* compute strings specific to individual A from population P<sub>A</sub> w.r.t. individual B from population P<sub>B</sub>

## Dependencies

C++11-compliant compiler (GCC 8.2 or newer), [ropebwt2](https://github.com/lh3/ropebwt2) and [htslib](https://github.com/samtools/htslib). For convenience, ropebwt2 and htslib are included in the repository.

## Download and Installation

```
git clone --recursive https://github.com/Parsoa/PingPong.git
cd PingPong 
cd ropebwt2 ; make ; cd ..
cd htslib ; make ; cd ..
make
```

You can now run PingPong by adding the clone directory to PATH. Because the package uses an internal clone of htslib, the shared objects will be in non-standard locations and have to be manually specified before running:

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/clone/dir/htslib
```

## How-To
Let's assume we have 3 samples A, B, and, C. To compute A-specific strings we have to:

1. Index samples B and C:
```
./PingPong index --binary --fastq /path/to/sample/B --index B.index.bin
./PingPong index --append B.index.bin --fastq /path/to/sample/C --index BC.index.fmd
```
2. Search for A-specific strings in the index
```
./PingPong search --index [B.index.bin] --fastq /path/to/sample/A --threads [nthreads]
```

The algorithm will output a [BED]() file named `subfreespecstrings.bed` with the list of A-specific strings. Each string is defined in terms of:
* identifier of the read it comes from
* starting position on the read
* length
* number of occurrences (we note that from this first pass, this number is always set to 1)

### PingPong Algorithm Usage
```
Usage: PingPong index [--binary] [--append /path/to/binary/index] --fastq /path/to/fastq --index /path/to/output/index

Optional arguments:
    -b, --binary          output index in binary format
    -a, --append          append to existing index (must be stored in binary)

Usage: PingPong search --index /path/to/index/file --fastq /path/to/fastq [--threads threads]

Optional arguments:
    --workdir             create output files in this directory (default:.)
    --overlap -1/0        run the exact algorithm (-1) or the relaxed one (0) (default:0)
    -t, --threads         number of threads (default:4)
```

##### Notes
* To append (`-a`) to an existing index, the existing index must be stored in binary format (`-b` option)
* An index built with `--binary` cannot be queried. Use `--binary` only for indices that are meant to be later appended to.
* The output file iscreated in the current directory (if `--workdir` is not set)
* Even when indexing a FASTA file, pass it with the `--fastq` option.

### Example

```
./PingPong index --binary --fastq example/father.fq --index example/father.fq.bin
./PingPong index --append example/father.fq.bin --fastq example/mother.fq --index example/index.fmd
./PingPong search --index example/index.fmd --fastq example/child.fq --overlap -1 --workdir example --threads 1
```

This will output strings that are specific to `child.fq` in `example/subfreespecstrings.bed`.

### Authors
For inquiries on this software please open an [issue](https://github.com/Parsoa/PingPong/issues) or contact either [Parsoa Khorsand](https://github.com/parsoa) or [Luca Denti](https://github.com/ldenti/).

### Citation

The PingPong is currently pending journal publication. The pre-print is availble below:

https://doi.org/10.1101/2021.03.23.436571
