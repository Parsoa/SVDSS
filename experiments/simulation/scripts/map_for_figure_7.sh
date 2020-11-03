#!/bin/bash
/share/hormozdiarilab/Codes/Stella/lib/bbmap/bbmap.sh threads=8 k=11 usequality=f idtag=t ref=alt.fa path=./bbmap in=short.fasta outm=short.sam
/share/hormozdiarilab/Codes/Stella/lib/bbmap/bbmap.sh threads=8 k=11 usequality=f idtag=t ref=../mother/alt.fa path=../mother/bbmap in=short.fasta outm=short.mother.sam
/share/hormozdiarilab/Codes/Stella/lib/bbmap/bbmap.sh threads=8 k=11 usequality=f idtag=t ref=../father/alt.fa path=../father/bbmap in=short.fasta outm=short.father.sam
cat short.sam | samtools sort -@ 8 -o short.sorted.bam
samtools index short.sorted.bam
cat short.mother.sam | samtools sort -@ 8 -o short.mother.sorted.bam
samtools index short.mother.sorted.bam
cat short.father.sam | samtools sort -@ 8 -o short.father.sorted.bam
samtools index short.father.sorted.bam
