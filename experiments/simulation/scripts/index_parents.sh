cd mother
cat chr*.fastq > reads.fastq
edit.sh
cd ..
cd father
cat chr*.fastq > reads.fastq
edit.sh
cd ..
cd mother-edited
/share/hormozdiarilab/Codes/Stella/stella pingpong index --binary --fastq reads_edited.fa > mother_edited.bwt.bin
/share/hormozdiarilab/Codes/Stella/stella pingpong index --append mother_edited.bwt.bin --fastq ../father-edited/reads_edited.fa > mother_edited+father_edited.bwt
