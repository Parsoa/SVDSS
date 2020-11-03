#!/bin/bash
#module load minimap2/2.17
echo $1, $2
RC=$(reverse_complement.sh $1)
echo $RC, "Reverse Complement"
grep $1 --before-context=1 --after-context=2 reads_edited.fa | sed '/^--$/d' >> $3/"$2".fa
grep $RC --before-context=1 --after-context=2 reads_edited.fa | sed '/^--$/d' >> $3/"$2".rc.fa
