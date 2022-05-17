#!/bin/sh
SVDSS_BIN=$1
fa=$2
bam=$3
wd=$4
threads=4

mkdir -p ${wd}

echo "################"
echo "### INDEXING ###"
echo "################"
${SVDSS_BIN} index --fastq ${fa} --index ${fa}.fmd

echo "#################"
echo "### SMOOTHING ###"
echo "#################"
${SVDSS_BIN} smooth --reference ${fa} --bam ${bam} --workdir ${wd} --threads ${threads}
samtools sort -T {output.bam}.sort-tmp ${wd}/smoothed.selective.bam > ${wd}/smoothed.selective.sorted.bam
samtools index ${wd}/smoothed.selective.sorted.bam

echo "#################"
echo "### SEARCHING ###"
echo "#################"
${SVDSS_BIN} search --index ${fa}.fmd --bam ${wd}/smoothed.selective.sorted.bam --threads ${threads} --workdir ${wd} --assemble

echo "###############"
echo "### CALLING ###"
echo "###############"
n=$(ls ${wd}/solution_batch_*.assembled.sfs | wc -l)
${SVDSS_BIN} call --reference ${fa} --bam ${wd}/smoothed.selective.sorted.bam --threads 4 --workdir ${wd} --batches ${n}
bcftools sort ${wd}/svs_poa.vcf > ${wd}/SVDSS.vcf