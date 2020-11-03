#!/bin/bash
P=$PWD
cd /share/hormozdiarilab/Codes/Stella
source venv3.6/bin/activate
cd $P
echo Plotting length/abundance correlation..
python ../scripts/plot.py corr solution_aggregated.fastq corr.png
echo Plotting SV length distribution..
python ../scripts/plot.py svlen de_novo.bed svlen.png
#echo Plotting short string abundances..
#python plot.py abund short_mapped.seqs.bed short_unmapped.bed abundance_short.png
#intersect long_mapped.bed ../child/de_novo.bed > long_mapped_sv.bed
#intersect long_mapped.bed ../child/de_novo.bed -v > long_mapped_non_sv.bed
#echo Plotting long string abundances..
#python plot.py abund long_mapped_sv.bed long_mapped_non_sv.bed abundance_long.png
#echo Plotting RepeatMasker ditribution for long strings..
#python plot.py repmask long_mapped.fasta long_mapped.fasta.out repeat_masker_long_mapped.png 0
#echo Plotting RepeatMasker ditribution for short strings..
#python plot.py repmask short_unmapped.fasta short_unmapped.fasta.out repeat_masker_short_unmapped.png 0
