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
    echo "[$(date)] Indexing.."
    ${SVDSS_BIN} index --reference ${fa} --index ${fa}.fmd
fi

if [ ! -f ${wd}/smoothed.bam ]
then
    echo "[$(date)] Smoothing.."
    ${SVDSS_BIN} smooth --reference ${fa} --bam ${bam} --threads ${threads} > ${wd}/smoothed.bam
    samtools index ${wd}/smoothed.bam
fi

if [ ! -f ${wd}/solution_batch_0.assembled.sfs ]
then
    echo "[$(date)] Computing SFSs.."
    ${SVDSS_BIN} search --index ${fa}.fmd --bam ${wd}/smoothed.bam --workdir ${wd} --threads ${threads}
fi

if [ ! -f ${wd}/SVDSS.vcf.gz ]
then
    echo "[$(date)] Calling SVs.."
    ${SVDSS_BIN} call --reference ${fa} --bam ${wd}/smoothed.bam --threads ${threads} --workdir ${wd} --min-cluster-weight ${w}
    bcftools sort ${wd}/svs_poa.vcf > ${wd}/SVDSS.vcf
fi
