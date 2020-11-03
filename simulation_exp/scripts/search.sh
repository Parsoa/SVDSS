/share/hormozdiarilab/Codes/Stella/src/cpp/stella/stella bedsearch --reference "$PWD/alt.fa" --bed $1
/share/hormozdiarilab/Codes/Stella/src/cpp/stella/stella bedsearch --reference "$PWD/../mother/alt.fa" --bed $1 | grep "Match" > found_mother.txt
/share/hormozdiarilab/Codes/Stella/src/cpp/stella/stella bedsearch --reference "$PWD/../father/alt.fa" --bed $1 | grep "Match" > found_father.txt
