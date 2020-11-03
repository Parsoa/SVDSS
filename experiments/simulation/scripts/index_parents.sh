cd mother
cat chr*.fastq > reads.fastq
edit.sh
cd ..
cd father
cat chr*.fastq > reads.fastq
edit.sh
cd ..
cd mother-edited
/share/hormozdiarilab/Codes/Stella/stella pingpong index -b reads_edited.fa > mother_edited.bwt.bin
/share/hormozdiarilab/Codes/Stella/stella pingpong index -a mother_edited.bwt.bin ../father-edited/reads_edited.fa > mother_edited+father_edited.bwt
