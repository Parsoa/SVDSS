#!/bin/sh

SVDSS_BIN=$1
fa=$2
bam=$3
wd=$4
threads=4
w=2

mkdir -p ${wd}

if [ ! -f ${fa}.fmd ]
then
    ${SVDSS_BIN} index --reference ${fa} --index ${fa}.fmd
fi

if [ ! -f ${wd}/smoothed.bam ]
then
    ${SVDSS_BIN} smooth --reference ${fa} --bam ${bam} --threads ${threads} > ${wd}/smoothed.bam
    samtools index ${wd}/smoothed.bam
fi

if [ ! -f ${wd}/specifics.txt ]
then
    ${SVDSS_BIN} search --index ${fa}.fmd --bam ${wd}/smoothed.bam --threads ${threads} > ${wd}/specifics.txt
fi

if [ ! -f ${wd}/calls.vcf ]
then
    ${SVDSS_BIN} call --reference ${fa} --bam ${wd}/smoothed.bam --sfs ${wd}/specifics.txt --threads ${threads} --min-cluster-weight ${w} | bcftools sort > ${wd}/calls.vcf
fi