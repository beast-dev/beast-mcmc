# big-3: BEAGLE against the MAT with the path and local blocks

Compares the Felsenstein likelihood (BEAGLE) with the MAT likelihood (integrated hits) using the path and
the local resampling blocks, on big-3 (899 BDBV genomes; HKY, strict clock, constant-size coalescent; the
XMLs of the beast-benchmark repository).

All three use the same data, model, priors and operators. The MAT patterns use `constantOnly="true"`
without `ambiguousConstant`, so that ambiguous states are treated exactly as BEAGLE treats them (17,194
variable sites rather than 2,405): with the ancestral states summed out, the MAT with integrated hits has
the same posterior as the Felsenstein likelihood, so the runs should agree once converged.

## Running

```
./make_xmls.py ../../../../beast-benchmark/Benchmarks/Benchmark-big xmls --chain-length 20000000 --log-every 5000
BEAST_CP="<classes>:<repo>/lib/*" ./run_benchmark.sh xmls runs "11 22" "mat-local mat-path baseline"
./analyse.py runs --burnin 0.25 --markdown results.md
```

The runs are made one after another, so that their timings are not disturbed by each other. BEAGLE uses
its CPU implementation with all 18 cores of the machine (an M5 Max); the MAT runs on one core. ESS per
minute is per minute of wall time after the burn-in.

## Results

(to come)
