# de-novo SVs from hifi reads

### Compilation
```
git clone https://github.com/Parsoa/Variation-Discovery.git
cd Variation-Discovery
 git submodule update --init --recursive ropebwt2
 cd ropebwt2
 make
 cd ..
 make
 ```

### How-To
```
# Create the index
./main index -b [father] > [father.index]
./main index -a [father.index] [mother] > [index]

# Search for child specific strings
./main sf3 [index] [child] [threads]
```

##### Note
* to append (`-a`) to an existing index, the existing index must be stored in binary format (`-b` option)
* to query the index, it must be stored in FMD format (default)
* the output is stored in the `solutions.fastq` file

### Example
```
./main index -b example/father.fq > example/father.fq.bin
./main index -a example/father.fq.bin example/mother.fq > example/index.fmd
./main sf3 example/index.fmd example/child.fq 1
 ```