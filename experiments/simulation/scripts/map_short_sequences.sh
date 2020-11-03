#!/bin/bash
module load minimap2/2.17
minimap2 -ax map-pb -t 8 $PWD/alt.fa short_mapped.fasta | samtools sort -@ 8 -o short_mapped.sorted.bam
minimap2 -ax map-pb -t 8 $PWD/alt.fa short_unmapped.fasta | samtools sort -@ 8 -o short_unmapped.sorted.bam
