# Experiments on real HiFi data

### Prerequisites
Install (everything is available on conda):
* bcftools
* bwa
* bbmap
* minimap2
* kmc
* ntedit
* bedtools
* snakemake
* biopython
* pysam
* numpy
* pandas
* matplotlib
* seaborn

Compile c++ scripts:
```
cd scripts
g++ -std=c++11 -O3 split.cpp -o split -lz
g++ -std=c++11 -O3 split_fa.cpp -o split_fa -lz
g++ -std=c++11 -I ~/miniconda3/include/ -o kmc2fa kmc2fa.cpp -L ~/miniconda3/lib -lkmc # or something similar
cd ..
```

### Data download
Get reference, contigs, and vcfs:
```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/hg38.no_alt.fa.gz
# extract chromosomes and save the new fasta (something like ref.fa)

wget -O hg.contigs1.fa http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20200509_HiFi_CHS-PUR-YRI_trio-hap-data/phased_assemblies/HG00733_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2_20200509.fasta
wget -O hg.contigs2.fa http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20200509_HiFi_CHS-PUR-YRI_trio-hap-data/phased_assemblies/HG00733_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2_20200509.fasta
cat hg.contigs1.fa hg.contigs2.fa | reformat.sh in=stdin.fa out=hg.contigs.fa uniquenames=t

wget -O na.contigs1.fa http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20200509_HiFi_CHS-PUR-YRI_trio-hap-data/phased_assemblies/NA19240_hgsvc_pbsq2-ccs_1000-pereg.h1-un.racon-p2.fasta
wget -O na.contigs2.fa http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20200509_HiFi_CHS-PUR-YRI_trio-hap-data/phased_assemblies/NA19240_hgsvc_pbsq2-ccs_1000-pereg.h2-un.racon-p2.fasta
cat na.contigs1.fa na.contigs2.fa | reformat.sh in=stdin.fa out=na.contigs.fa uniquenames=t

# Note: we used freeze 4 - I'll update the link once the vcfs will be available
wget -O snps.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20200814_Freeze3_merged_PAV_callset/freeze3.snv.alt.vcf.gz
wget -O indels.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20200814_Freeze3_merged_PAV_callset/freeze3.indel.alt.vcf.gz
wget -O svs.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20200814_Freeze3_merged_PAV_callset/freeze3.sv.alt.vcf.gz
```

### Samples
Get the samples
* [HG00733](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20190925_PUR_PacBio_HiFi/) (all HG00733 runs except the LAMBDA contaminated one)
* [NA19240](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20191005_YRI_PacBio_NA19240_HiFi/)
concatenate them (e.g., `zcat * > `) and correct them with `ntEdit`:
```
/usr/bin/time -v nthits -t {threads} -k 31 -c 2 -b 36 --outbloom -p {outdir}/ {sample}
/usr/bin/time -v ntedit -t {threads} -f {sample} -r {outdir}/_k31.bf -b {outdir}/
```

### Merge and clean VCFs
```
tabix -p vcf snps.vcf.gz
tabix -p vcf indels.vcf.gz
tabix -p vcf svs.vcf.gz

# Concatenate VCFs and remove records with missing genotypes
bcftools concat -a snps.vcf.gz indels.vcf.gz svs.vcf.gz > all.vcf
python3 scripts/remove_missing_gt.py all.vcf > all.known.vcf
bgzip all.known.vcf
tabix -p vcf all.known.vcf.gz

# Try to create the consensus: it will find overlapping variations
bcftools consensus -s HG00733 -H 1 -f ref.fa all.known.vcf.gz > /dev/null 2> log.1
bcftools consensus -s HG00733 -H 2 -f ref.fa all.known.vcf.gz > /dev/null 2> log.2

# Create the new VCF removing overlapping variations
python3 scripts/filter_removed.py all.known.vcf.gz log.1 log.2 > all.known.nonoverlapping.vcf
bgzip all.known.nonoverlapping.vcf
tabix -p vcf all.known.nonoverlapping.vcf.gz
```

### Run Experiments
1. set variables in the `config.yml` file (notes: everything will be store in the `outdir` folder. `codedir` is the local path to this repo)
2. run snakemake
