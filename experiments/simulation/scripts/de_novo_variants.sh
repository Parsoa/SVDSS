intersect present.bed ../mother/present.bed -v > de_novo_mother.bed | wc -l
intersect de_novo_mother.bed ../father/present.bed -v > de_novo.bed | wc -l
