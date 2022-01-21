![C/C++ CI](https://github.com/Parsoa/PingPong/workflows/C/C++%20CI/badge.svg)

# SVDSS: Structural Variant Discovery from Sample-specific Strings

SVDSS is a novel method for discovery of structural variants in accurate long reads (e.g PacBio HiFi) using sample-specific strings (SFS)

Efficient computation of A-specific string w.;.t. a set {B,C,...,Z} of other long reads samples. A A-specific string is a string which occurs only in sample A and not in the others. 

**Note:** if you are looking for the _biorxiv_ implementation, please refer to [this release](https://github.com/Parsoa/PingPong/releases/tag/v1.0.0-pingpong).

##### Use-Cases

* compute strings specific to child w.r.t. parents
* compute strings specific to individual A from population P<sub>A</sub> w.r.t. individual B from population P<sub>B</sub>

## Dependencies

C++11-compliant compiler (GCC 8.2 or newer). The following libraries are need to build and run SVDSS, all are included as submodules:

* [htslib](https://github.com/samtools/htslib) built with [libdeflate](https://github.com/ebiggers/libdeflate) for BAM processing.
* [ksw2](https://github.com/lh3/ksw2) for FASTA and FASTQ processing.
* [ropebwt2](https://github.com/lh3/ropebwt2) for FMD index creation and querying.
* [abPOA](https://github.com/yangao07/abPOA) for POA computation.
* [parasail](https://github.com/jeffdaily/parasail) for local alignment of POA consensus.
* [rapidfuzz] (https://github.com/maxbachmann/rapidfuzz-cpp) for string similarity computation.
* [interval-tree] (https://github.com/5cript/interval-tree) for variant overlap detection and clustering.

All libraries are included as submodules in the repository. Note that SVDSS requires `htslib` built with `libdeflate` for performance reasons, so you still need to build the local version even if you have it installed globally. We refer to the Makefile for additional details about building dependencies.

## Download and Installation

You need to clone the repostiory with `--recursive` for all the dependencies to be downloaded. Each dependency has to be built separately:

```
git clone --recursive https://github.com/Parsoa/PingPong.git
cd SVDSS 
git submodule update --init --recursive

cd parasail
mkdir build ; cd build
cmake ..
make
cd ../..

cd rapidfuzz-cpp
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .

cd abPOA
make
cd ..

cd libdeflate
make
cd ..

cd htslib
./configure --with-libdeflate
make
cd ..

make
```

You need add some of the paths to your envinronment for the libraries to link properly at runtime:

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/cloned/repo/htslib:/path/to/cloned/repo/parasail/build:/path/to/cloned/repo/libdeflate
```

## General Usage

```
    Index genome:
        SVDSS index [--binary] [--append /path/to/binary/index] --fastq /path/to/fastq/file --index /path/to/output/index/file
        Optional arguments: 
            -b, --binary                output index in binary format
            -a, --append                append to existing index (must be stored in binary). DON'T pass this option for building an index you want to use directly.
    Extract SFS from BAM/FASTQ/FASTA files:
        SVDSS search --index /path/to/index --fastq/--bam /path/to/input --workdir /output/directory
        Optional arguments: 
            --assemble                automatically runs SVDSS assemble on output
    Assmble SFS into superstrings
        SVDSS assemble --workdir /path/to/.sfs/files --batches /number/of/SFS/batches
    Reconstruct sample
        SVDSS reconstruct --workdir /output/file/direcotry --bam /path/to/input/bam/file --reference /path/to/reference/genome/fasta
        Optional arguments: 
            --selective                ignores reads with high mismatch rate
    Call SVs:
        SVDSS call --workdir /path/to/assembled/.sfs/files --bam /path/to/input/bam/file --reference /path/to/reference/genome/fasta
        Optional arguments: 
            --clipped                calls SVs from clipped SFS only.
            --min-cluster-weight                minimum number of supporting superstrings for a call to be reported.
            --min-sv-length                minimum length of reported SVs. Default is 25. Values < 25 are ignored.
    General options: 
        --threads                sets number of threads, default 4.
```

## Usage Guide

SVDSS requires as input the BAM file of the sample to be genotyped, a reference genome in FASTA format. To genotype a sample we need to perform the following steps:

1. Build FMD index of reference genome (`SVDSS index`)
2. Reconstruct the input BAM file (`SVDSS reconstruct`)
3. Extract SFS from reconstructed BAM file (`SVDSS search`)
4. Assemble SFS into superstrings (`SVDSS assemble`)
5. Genotype SVs from the assembled superstrings (`SVDSS call`)

Figure below shows the full pipeline of commands that needs to be run:

![SVDSS's pipeline](docs/Pipeline.png)

Note that you can reuse the index from step 1 for any number of samples genotyped against the same reference genome.

In the guide below we assume we are using the reference genome file `GRCh38.fa` and the input BAM file `sample.bam`. We assume both files are present in the working directory. All of SVDSS steps must be run in the same directory so we always pass `--workdir $PWD` for every command.

### Index reference genome

The FMD index is the same as from PingPong:

```
SVDSS index --fastq GRCh38.fa --index GRCh38.bwt
```

The `--index` option specifies the output file name.

### Reconstruct target sample

Reconstruction removes  nearly all SNPs, small indels and sequencing errors from reads. This results in smaller number of SFS being extracted and increases the relevance of extracted SFS to SV discovery significantly. To reconstruct the sample run: 

```
SVDSS reconstruct --bam sample.bam --workdir $PWD --reference GRCh38.fa --threads 16
```

This produces a file named `reconstructed.selective.bam`. This file is sorted in the same order as the input file, however it needs to be indexed again with `samtools index`. The command also produces two files `recontructed_reads.txt` and `ignored_reads.txt` in `workdir` that contains the ids of reads that were reconstructed and ids of reads that didn't have any large (> 20bp) indels in their alignemnts. This information is used by the next step.

### Extract SFS from target sample

To extract SFS run:

```
SVDSS search --index GRCh38.bwt --bam reconstructed.selective.bam --workdir $PWD
```

This step produces a number of `solution_batch_<i>.sfs` files. These files include the coordinates of SFS relative to the reads they were extracted from.

### Assemble SFS into superstrings

To reduce redundancy, overlapping SFS on each reads are merged. Simply run:

```
SVDSS assemble --workdir $PWD --batches N
```

Here `N` is the number of files produces by the previous step. The output will be the same number of `solution_batch_<i>.assembled.sfs`.

### Call SVs

We are now ready to call SVs. Run:

```
SVDSS call --reference GRCh38.fasta --bam reconstructed.selective.bam --workdir $PWD --batches N
```

You can filter the reported SVs by passing the `--min-sv-length` and `--min-cluster-weight` options. These options control the minimum length and minimum number of supporting superstrings for the reported SVs. Higher values for `--min-cluster-weight` will increase precision at the cost of reducing recall. For a 30x coverage sample, `--min-cluster-weight 4` produced the best results in our experiments.

This commands output two files: `svs_poa.vcf` that includes the SV calls and `poa.sam` which includes alignments of POA contigs to the reference genome.

### Authors
For inquiries on this software please open an [issue](https://github.com/Parsoa/SVDSS/issues) or contact either [Parsoa Khorsand](https://github.com/parsoa) or [Luca Denti](https://github.com/ldenti/).

### Citation

Our manuscript describing the SFS idea and the PingPong algorithm titled "Comparative genome analysis using sample-specific string detection in accurate long reads" is now published in Bioinformatics Advances:

https://doi.org/10.1093/bioadv/vbab005
