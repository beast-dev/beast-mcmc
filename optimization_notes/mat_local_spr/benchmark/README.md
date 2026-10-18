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
minute is per minute of wall time after the burn-in. The tree-related parameters burn in slowly from the
starting tree: use a burn-in of 50% for runs of this length.

## Results

Two seeds per method, 20 million states each (October 2026, `operator_experiment` at b7ec573e5), on an
M5 Max. The full tables are in [results.md](results.md) (50% burn-in).

| Method | Time per million states | Cores | ESS per minute: root height | tree length | population size | clock rate | coalescent | kappa | lowest in basket |
|---|---|---|---|---|---|---|---|---|---|
| BEAGLE (baseline) | 2.9 min | 18 | 2.1 | 1.6 | 2.6 | 2.5 | 2.0 | 36 | 1.6 |
| MAT, path block | 2.8 min | 1 | 2.0 | 0.8 | 1.1 | 1.2 | 0.7 | 41 | 0.7 |
| MAT, local block | 1.15 min | 1 | 7.7 | 3.3 | 4.6 | 5.2 | 3.6 | 94 | 3.3 |

- **Agreement.** With 50% burn-in all three methods agree: every difference from the baseline is within
  1.5 Monte Carlo standard errors, as are the differences between the replicates of each method and between
  the halves of each run. (With 25% burn-in both MAT methods still showed slightly longer trees, larger
  populations and lower clock rates than the baseline, by 2-3.5 standard errors: residual burn-in, which
  disappears with the longer burn-in.) The MAT with integrated hits therefore samples the same posterior as
  the Felsenstein likelihood, as it should.
- **Speed.** The local block runs 2.5 times as many states per minute as BEAGLE (on one core rather than
  18) and 2.4 times as many as the path block.
- **Mixing per state.** The local block mixes about as well per state as BEAGLE on the slow parameters
  (tree length, population size, clock rate), despite its smaller resampling block, so its ESS per minute is
  about twice to three times BEAGLE's: 3.3 against 1.6 for the slowest parameter. The path block is no
  faster than BEAGLE and mixes less well per state on these parameters, so it has about half BEAGLE's ESS
  per minute.
- **Caveats.** ESS values of 10-150 (the tree-related parameters) are themselves uncertain, and they
  differ between replicates by up to a factor of three; the comparison is between two runs per method.
  BEAGLE uses all 18 cores, the MAT one: per core, the local block is about 40 times as efficient.
