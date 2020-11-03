#!/bin/bash
module load minimap2/2.17
minimap2 -ax map-pb -t 8 ../child/alt.fa long_not_covering.fasta | samtools sort -@ 8 -o long_not_covering.bam
samtools index long_not_covering.bam
/share/hormozdiarilab/Codes/Stella/lib/bbmap/bbmap.sh threads=8 k=11 usequality=f idtag=t ref=../child/alt.fa path=./bbmap in=short_mapped.fasta outm=short_mapped.sam
cat short_mapped.sam | samtools sort -@ 8 -o short_mapped.sorted.bam
/share/hormozdiarilab/Codes/Stella/lib/bbmap/bbmap.sh threads=8 k=11 usequality=f idtag=t ref=../child/alt.fa path=./bbmap in=short_unmapped.fasta outm=short_unmapped.sam
cat short_unmapped.sam | samtools sort -@ 8 -o short_unmapped.sorted.bam
