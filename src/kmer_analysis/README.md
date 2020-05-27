# Approach 1.1
1. index the parents' samples with an FM-index
2. count the kmers on the child sample 
3. search each kmer in the index counting how many kmers are child-specific

### Compilation
```
git submodule update --init --recursive sdsl-lite KMC
cd sdsl-lite
bash install.sh ..
cd ..
cd KMC
make
cd ..
make
```

##### Examples
```
# Step 1: create the index (1 thread)
./main index example/ERR3840017.fq.gz example/ERR3861388.fq.gz example/index.fm9

# Step 2: run KMC (4 threads, -t)
mkdir -p kmc_tmp
./KMC/bin/kmc -b -k17 -t4 example/ERR3861386.fq.gz example/ERR3861386.kmc kmc_tmp
# rm -r kmc_tmp

# Step 3: count child-specific kmers (4 threads, 3rd argument)
./main search example/index.fm9 example/ERR3861386.kmc 4 > example/cskmers.txt
```

### Analysis
1. set data folder, samples, chromosomes, and kmer sizes in `config.yaml`
2. run snakemake to:
    1. download the data (reference and samples)
    2. align the samples
    3. extract reads from specific chromosomes
    4. run approach1.1


##### FIXME
* create a conda environment with all the dependencies (minimap2, samtools...)
