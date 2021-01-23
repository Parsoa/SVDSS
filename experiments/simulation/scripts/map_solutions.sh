#!/bin/bash
module load minimap2/2.17
module load samtools
cat solution_aggregated.fastq | awk 'BEGIN {counter=0; header=""} {if (NR % 4 == 1) {header = $0} else if (NR % 4 == 2) {if (length($0) < 500) {counter = 1} else {counter = 0; print header; print $0} } else { if (counter == 0) {print $0} } }' > solution_aggregated.long.fastq
cat solution_aggregated.fastq | awk 'BEGIN {counter=0; header=""} {if (NR % 4 == 1) {header = $0} else if (NR % 4 == 2) {if (length($0) > 500) {counter = 1} else {counter = 0; print header; print $0} } else { if (counter == 0) {print $0} } }' > solution_aggregated.short.fastq
minimap2 -ax map-pb -t 8 $PWD/alt.fa solution_aggregated.long.fastq | samtools sort -@ 8 -o solution.aggregated.long.bam
samtools index solution.aggregated.long.bam
/share/hormozdiarilab/Codes/Stella/lib/bbmap/bbmap.sh threads=8 k=11 ambiguous=all usequality=f idtag=t ref=../child/alt.fa path=./bbmap in=solution_aggregated.short.fastq outm=solution.aggregated.short.sam
cat solution.aggregated.short.sam | samtools sort -@ 8 -o solution.aggregated.short.bam
