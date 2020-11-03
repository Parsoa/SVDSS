export PATH=$PATH:/share/hormozdiarilab/Codes/Stella/lib/ntEdit
export PATH=$PATH:/share/hormozdiarilab/Codes/Stella/lib/ntHits
P=$PWD
cd ..
mkdir "$P"-edited
cd "$P"-edited
nthits -t 8 -k 31 -c 2 -b 36 --outbloom -p reads "$P"/reads.fastq
ntedit -t 8 -m 2 -i 4 -d 5 -X 0.5 -Y 0.5 -f "$P"/reads.fastq -r reads_k31.bf -b reads
convert_fasta_to_fastq.py
cd ..
