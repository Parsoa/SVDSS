
Get data:
```
mkdir {root}
cd {root}

# callset
wget <vcf>
wget <vcf>.tbi

# reference
wget <ref>
samtools faidx hg38.no_alt.fa.gz chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > hg38.chroms.fa
samtools faidx hg38.chroms.fa

# assemblies
wget
wget

wget
wget

# samples
wget + cat
wget + cat
```

Setup conda:
```
conda install snakemake-minimal bedtools bcftools samtools biopython bbmap ntedit pysam seaborn matplotlib
```

Compile scripts:
```
g++ -std=c++11 -O3 -o aggregate aggregate.cpp -lz
g++ -std=c++11 -O3 -o filter filter.cpp -lz
```

RepeatMasker:
```
python3 scripts/plot.py repmask <fa> <rm.out> <png>
```
