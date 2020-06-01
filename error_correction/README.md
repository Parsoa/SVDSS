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
./cutsample [sample] [kmc_prefix] [threshold]
```
