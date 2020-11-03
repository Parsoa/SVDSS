mkdir short
samtools view solution.aggregated.short.bam | awk '{print $1, $3, $4, $10}' | awk '{if ($3 != 0) {print $0} }' > ./short/solution.aggregated.short.mapped.txt
samtools view solution.aggregated.short.bam | awk '{print $1, $3, $4, $10}' | awk '{if ($3 == 0) {print $0} }' > ./short/solution.aggregated.short.unmapped.txt
cd short 
split -d -n l/1000 --additional-suffix=.txt solution.aggregated.short.mapped.txt solution.mapped.batch_
split -d -n l/1000 --additional-suffix=.txt solution.aggregated.short.unmapped.txt solution.unmapped.batch_
cd ..
mkdir long
samtools view solution.aggregated.long.bam | awk '{print $1, $3, $4, $10}' | awk '{if ($3 != 0) {print $0} }' > ./long/solution.aggregated.long.mapped.txt
samtools view solution.aggregated.long.bam | awk '{print $1, $3, $4, $10}' | awk '{if ($3 == 0) {print $0} }' > ./long/solution.aggregated.long.unmapped.txt
cd long 
split -d -n l/1000 --additional-suffix=.txt solution.aggregated.long.mapped.txt solution.mapped.batch_
split -d -n l/1000 --additional-suffix=.txt solution.aggregated.long.unmapped.txt solution.unmapped.batch_
