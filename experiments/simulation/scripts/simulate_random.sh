P=$PWD
echo $P
stella.sh simulate --bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/HGSV/Unified/HG00733.unified.all.bed /share/hormozdiarilab/Codes/NebulousSerendipity/data/HGSV/Inversions/HG00733.INV.bed --workdir "$PWD" --reference /share/hormozdiarilab/Data/ReferenceGenomes/Hg38/GRC38.fasta --random --coverage $1
cd $P
index_parents.sh
cd "$P"/child
index_reference.sh
cat chr*.fastq > reads.fastq
edit.sh
cd "$P"/child-edited
solve.sh
map_solutions.sh
parse_solution.sh
find.sh
