# BEAGLE against the MAT with the path and local blocks

Compares the Felsenstein likelihood (BEAGLE) with the MAT likelihood (integrated hits) using the path and
the local resampling blocks, with and without the Gibbs operator on the mutational history
(`mutationHistoryGibbsOperator`), on two datasets from the beast-benchmark repository:

- **big-3:** 899 BDBV genomes; HKY, strict clock, constant-size coalescent. The speed benchmark, but its
  tree parameters mix slowly (ESS 13-150 in 20 million states).
- **small-2:** 165 influenza A HA sequences; JC69, strict clock, constant-size coalescent. The validation:
  its tree parameters mix well (ESS 300-1900 in 50 million states), so a small bias would show.

All methods use the same data, model, priors and operators. The MAT patterns use `constantOnly="true"`
without `ambiguousConstant`, so that ambiguous states are treated exactly as BEAGLE treats them (on big-3,
17,194 variable sites rather than 2,405): with the ancestral states summed out, the MAT with integrated
hits has the same posterior as the Felsenstein likelihood, so the runs should agree once converged.

## Running

```
./make_xmls.py ../../../../beast-benchmark/Benchmarks/Benchmark-big xmls --chain-length 20000000 --log-every 5000
BEAST_CP="<classes>:<repo>/lib/*" ./run_benchmark.sh xmls runs "11 22 33 44" "mat-local-gibbs mat-local baseline mat-path-gibbs"
./analyse.py runs --burnin 0.5 --mat-reference mat-path-gibbs --markdown results.md

./make_validation_xmls.py ../../../../beast-benchmark/Benchmarks/Benchmark-small validation-xmls --chain-length 50000000
BEAST_CP="<classes>:<repo>/lib/*" ./run_benchmark.sh validation-xmls validation-runs "11 22 33 44" "mat-local-gibbs mat-local baseline"
./analyse.py validation-runs --burnin 0.1 --mat-reference mat-local-gibbs --markdown validation-results.md
```

The runs are made one after another, so that their timings are not disturbed by each other. BEAGLE uses
its CPU implementation with all 18 cores of the machine (an M5 Max); the MAT runs on one core. ESS per
minute is per minute of wall time after the burn-in. Keep the machine awake (`caffeinate -s` on mains
power): wall times that include sleep are useless.

Besides the basket of model parameters, `analyse.py` reports for the MAT runs their own likelihood, the
number of mutations in the sampled history (`mat.mutations`) and its **excess** over the Fitch parsimony
score of the sampled tree (from the .trees file, `tree_parsimony.py`). A history that mixes slowly keeps
more mutations than its tree needs: an excess that is high and differs between replicates. The tree
parameters alone did not show this on big-3.

## Results

October 2026, `operator_experiment` at 7d0e64742 (the BEAGLE runs of big-3 seeds 11 and 22 at b7ec573e5,
which differs only in the MAT), four seeds per method. Full tables: [results.md](results.md) (big-3, 50%
burn-in) and [validation-results.md](validation-results.md) (small-2, 10% burn-in).

### The local block needs the Gibbs operator

On small-2, posterior means of the four seeds:

| Method | Tree length | Root height | Population size | Clock rate (×10⁻³) | Excess mutations |
|---|---|---|---|---|---|
| BEAGLE | 62.48-62.68 | 6.20-6.21 | 3.70-3.73 | 3.72-3.74 | |
| MAT local + Gibbs | 62.16-62.69 | 6.18-6.22 | 3.69-3.73 | 3.72-3.75 | 0.31-0.34 |
| MAT local, no Gibbs | 62.75-67.63 | 6.21-6.45 | 3.73-4.17 | 3.50-3.72 | 0.85-4.84 |

- **With the Gibbs operator the MAT samples the same posterior as BEAGLE:** pooled over seeds, every
  parameter is within 1.1 Monte Carlo standard errors of BEAGLE, and the seeds agree with each other.
- **Without it the local block does not.** The local block resamples a node's states only when a subtree
  move regrafts next to it, so excess mutations elsewhere in the history are almost never removed: each
  run keeps its own excess for all 50 million states, and the tree parameters follow it (more mutations
  fit longer trees, larger populations and a lower clock rate). Pooled, tree length is 16 standard errors
  from the Gibbs runs, which match BEAGLE. The local block is only correct together with the Gibbs operator.
- The same happens on big-3 (excess 0.2-5.1 without the operator, 0.04-0.14 with it), but there the tree
  parameters mix too slowly for it to show: the earlier agreement of the plain local block with BEAGLE
  (two seeds) was luck.

### Speed

| | Method | Time per million states | Cores | ESS per minute: root height | tree length | population size | clock rate | coalescent |
|---|---|---|---|---|---|---|---|---|
| small-2 | BEAGLE | 0.70 min | 18 | 55 | 28 | 31 | 35 | 24 |
| | MAT local + Gibbs | 0.30 min | 1 | 86 | 43 | 58 | 56 | 44 |
| big-3 | BEAGLE | 3.0 min | 18 | 3.5 | 1.4 | 2.1 | 2.3 | 1.6 |
| | MAT local + Gibbs | 0.78 min | 1 | 10.1 | 5.8 | 8.6 | 9.5 | 6.7 |
| | MAT path + Gibbs | 2.4 min | 1 | 2.7 | 1.5 | 2.2 | 2.3 | 1.7 |

- The Gibbs operator (weight about 10% of the SubtreeLeap weight) costs little: on big-3 the local block
  runs at 0.78 min per million states with it and 0.76 without.
- With it, the local block has about 1.5-2 times BEAGLE's ESS per minute on small-2 and about 3-4 times on
  big-3, on one core against 18.
- **Caveats.** On small-2 the ESS per minute is from seeds 11 and 22 only: the machine slept during the
  later seeds, or ran other work, so their wall times are not usable (their samples are). On big-3 the
  wall times of the path block seeds 33 and 44 include sleep and are replaced by the mean of seeds 11 and
  22 (see runs/BUILD.txt). On big-3 the tree parameters have ESS of 20-150 per run, and the four seeds of
  each MAT method with the Gibbs operator differ by up to 5.6 standard errors in tree length (pooled, about
  2% longer trees than BEAGLE). The traces show seeds held in different regions of tree space for millions
  of states, and small-2, where the trees mix well, shows no bias: slow tree mixing, but big-3 needs much
  longer runs for a strict comparison.
