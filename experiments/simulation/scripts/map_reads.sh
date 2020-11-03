#!/bin/bash
module load minimap2/2.17
minimap2 -ax map-pb -t 8 /share/hormozdiarilab/Data/ReferenceGenomes/Hg38/GRC38.fasta reads.fastq > reads.sam
cat reads.sam | samtools sort -@ 8 -o reads.sorted.bam
samtools index reads.sorted.bam
