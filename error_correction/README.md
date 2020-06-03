# Error "Correction"

```
cd KMC
make
cd ..
make
#mkdir -p kmc_tmp
#./KMC/bin/kmc -k3 test/sample.fq test/sample.fq.k3 kmc_tmp
./cutsample test/sample.fq test/sample.k3 5
```

```
./cutsample [sample] [kmc_prefix] [threshold] > [fragmented_sample]
```

```
python3 print_len_fastx.py [sample] > [sample.len]
python3 plot_fraglen_dist.py [real_sample.len] [fragment.len] [out_prefix]
```

```
python3 print_error_len.py [fragmented_sample] > [fragmented_sample.gaps]
python3 plot_gaplen_dist.py [fragmented_sample.gaps] # removes SNPs by default
```
