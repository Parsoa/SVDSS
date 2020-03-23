# Approach 1.1
1. index the parents' samples with an FM-index
2. count the kmers on the child sample 
3. search each kmer in the index counting how many kmers are "unique" to the child

### Compilation
```
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
mkdir -p kmc_tmp
./KMC/bin/kmc -k17 example/ERR3861386.fq.gz example/ERR3861386.kmc kmc_tmp
./main example/ERR3840017.fq.gz example/ERR3861388.fq.gz example/ERR3861386.kmc
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
