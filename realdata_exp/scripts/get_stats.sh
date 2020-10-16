#!/bin/bash

# TODO remove duplicate code

wdir=$1

pp=$1/pp
step1=$1/step1
step2=$1/step2

### ping-pong ###
echo "#specific: $(cut -f 5 -d' ' ${pp}/aggregate.log)"
echo "#specific (ab>=2): $(sed -n '1~4p' ${pp}/sspecific.fq | wc -l)"
echo "#specific (ab>=5): $(sed -n '1~4p' ${pp}/sspecific.fq | grep -v -P '#2$|#3$|#4' | wc -l)"
echo "#specific (ab>=5; len<=500): $(sed -n '1~4p' ${pp}/sspecific.minab_5.maxl_500.fq | wc -l)"

echo ""
echo ""

### step1 ###
HG=${step1}/sspecific.minab_5.maxl_500.cleared.HG00733.sam
NA=${step1}/sspecific.minab_5.maxl_500.cleared.NA19240.sam

echo "Tot alignments to HG contigs: $(samtools view -c ${HG})"
echo "Primary alignments to HG contigs: $(samtools view -F 256 -c ${HG})"
echo "Aligned perfectly to HG contigs: $(python3 extract_al_by_err.py ${HG} 0 none |& cut -f 4 -d' ')"
echo "Aligned with 1 err to HG contigs: $(python3 extract_al_by_err.py ${HG} 1 none |& cut -f 4 -d' ')"

echo ""

echo "Tot alignments to NA contigs: $(samtools view -c ${NA})"
echo "Primary alignments to NA contigs: $(samtools view -F 256 -c ${NA})"
echo "Aligned perfectly to NA contigs: $(python3 extract_al_by_err.py ${NA} 0 none |& cut -f 4 -d' ')"
echo "Aligned with 1 err to NA contigs: $(python3 extract_al_by_err.py ${NA} 1 none |& cut -f 4 -d' ')"

echo ""
echo ""

HGhap=${step2}/sspecific.minab_5.maxl_500.cleared.HG00733.sam
echo "Tot alignments to HG haps: $(samtools view -c ${HGhap})"
echo "Primary alignments to HG haps: $(samtools view -F 256 -c ${HGhap})"
echo "Aligned perfectly to HG haps: $(python3 extract_al_by_err.py ${HGhap} 0 none |& cut -f 4 -d' ')"
echo "Aligned with 1 err to HG haps: $(python3 extract_al_by_err.py ${HGhap} 1 none |& cut -f 4 -d' ')"

echo ""
echo ""

HGunc=${step2}/uncovering/sspecific.HG00733.sam
NAunc=${step2}/uncovering/sspecific.NA19240.sam

echo "Tot uncovering to HG contigs: $(samtools view -c ${HGunc})"
echo "Primary uncovering to HG contigs: $(samtools view -F 256 -c ${HGunc})"
echo "Uncovering aligned perfectly to HG contigs: $(python3 extract_al_by_err.py ${HGunc} 0 none |& cut -f 4 -d' ')"
echo "UNcovering aligned with 1 err to HG contigs: $(python3 extract_al_by_err.py ${HGunc} 1 none |& cut -f 4 -d' ')"

echo ""

echo "Tot uncovering to NA contigs: $(samtools view -c ${NAunc})"
echo "Primary uncovering to HG contigs: $(samtools view -F 256 -c ${NAunc})"
echo "Uncovering aligned perfectly to HG contigs: $(python3 extract_al_by_err.py ${NAunc} 0 none |& cut -f 4 -d' ')"
echo "UNcovering aligned with 1 err to HG contigs: $(python3 extract_al_by_err.py ${NAunc} 1 none |& cut -f 4 -d' ')"
