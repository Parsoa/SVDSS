#!/bin/bash
module load minimap2/2.17
/share/hormozdiarilab/Codes/Stella/main sf3 ../mother/mother_edited+father_edited.bwt reads.fastq 48
cat solution_batch_* > solution.fastq
/share/hormozdiarilab/Codes/Stella/src/cpp/stella/stella aggregate --workdir $PWD --threads 48
