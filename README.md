![C/C++ CI](https://github.com/Parsoa/PingPong/workflows/C/C++%20CI/badge.svg)

# Sample-specific string detection from accurate long reads
Efficient computation of A-specific string w.r.t. a set {B,C,...,Z} of other long reads samples. A A-specific string is a string which occur only in sample A and not in the others. 

##### Use-Cases
* compute strings specific to child w.r.t. parents
* compute strings specific to individual A from population P<sub>A</sub> w.r.t. individual B from population P<sub>B</sub>

## Dependencies

C++11-compliant compiler, [ropebwt2](https://github.com/lh3/ropebwt2) and [htslib](https://github.com/samtools/htslib). For convenience, ropebwt2 and htslib are included in the repository.

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

The algorithm will output a set of files named `solution_batch_<number>.fastq` with the list of A-specific strings. This is done to keep the memory usage constant. A string maybe repeated in multiple files and abundances from all files should be aggregated to represent final numbers.

3. Aggregate output to a single fastq file:
```
./PingPong aggregate --workdir /path/to/string/batches --threads <threads> --cutoff <minimum abundance for strings> --batches <number of batches>
```

The `--batches` option should be set to the number of `solution_batch_<number>.fastq` files produced, e.g if the last file is `solution_batch_5.fastq` then `--batches` should be set to 6 (batch file index starts from zero).

You can also directly aggregate output from the search step by adding `--aggregate --cutoff <value>` to the command in step 2 and then skip step 3.

### PingPong Algorithm Usage
```
Usage: PingPong index [--binary] [--append index] --fastq /path/to/fastq/file --index /path/to/output/index/file [--threads threads]

Optional arguments:
    -b, --binary          output index in binary format
    -a, --append          append to existing index (must be stored in binary)
    -t, --threads         number of threads (default:1)

Usage: PingPong search [--index /path/to/index/file] [--fastq fastq] [--threads threads]

Optional arguments:
    --aggregate         aggregate ouputs directly.
    --cutof             sets cutoff for minimum string abundance (tau)
```

##### Notes
* To append (`-a`) to an existing index, the existing index must be stored in binary format (`-b` option)
* An index built with `--binary` cannot be queried. Use `--binary` only for indices that are meant to be later appended to.
* The search output is stored in multiple `solution_batch_*.fastq` files (created in the current directory)
* Even when indexing a FASTA file, pass it with the `--fastq` option.

### Example

```
./PingPong index --binary --fastq example/father.fq --index example/father.fq.bin
./PingPong index --append example/father.fq.bin --fastq example/mother.fq --index example/index.fmd
./PingPong search --index example/index.fmd --fastq example/child.fq --threads 4
 ```
 
 ### Authors
 For inquiries on this software please open an [issue](https://github.com/Parsoa/PingPong/issues) or contact either [Parsoa Khorsand](https://github.com/parsoa) or [Luca Denti](https://github.com/ldenti/).
