#!/bin/bash
module load minimap2/2.17
minimap2 -ax map-pb -t 8 ../child/alt.fa long.fasta | samtools sort -@ 8 -o long.bam
samtools index long.bam
minimap2 -ax map-pb -t 8 ../mother/alt.fa long.fasta | samtools sort -@ 8 -o long.mother.bam
samtools index long.mother.bam
minimap2 -ax map-pb -t 8 ../father/alt.fa long.fasta | samtools sort -@ 8 -o long.father.bam
samtools index long.father.bam

/share/hormozdiarilab/Codes/Stella/lib/bbmap/bbmap.sh threads=8 k=11 usequality=f idtag=t ref=../child/alt.fa path=./bbmap in=short.fasta outm=short.sam
/share/hormozdiarilab/Codes/Stella/lib/bbmap/bbmap.sh threads=8 k=11 usequality=f idtag=t ref=../mother/alt.fa path=../mother/bbmap in=short.fasta outm=short.mother.sam
/share/hormozdiarilab/Codes/Stella/lib/bbmap/bbmap.sh threads=8 k=11 usequality=f idtag=t ref=../father/alt.fa path=../father/bbmap in=short.fasta outm=short.father.sam
cat short.sam | samtools sort -@ 8 -o short.sorted.bam
cat short.mother.sam | samtools sort -@ 8 -o short.mother.sorted.bam
cat short.father.sam | samtools sort -@ 8 -o short.father.sorted.bam
