echo $@
/share/hormozdiarilab/Codes/Stella/src/cpp/stella/stella find --workdir $PWD --reference $PWD/alt.fa --threads 48 --bed $PWD/present.bed --type long  "$@"
/share/hormozdiarilab/Codes/Stella/src/cpp/stella/stella find --workdir $PWD --reference $PWD/alt.fa --threads 48 --bed $PWD/present.bed --type short "$@"
