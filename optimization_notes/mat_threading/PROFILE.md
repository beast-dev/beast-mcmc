# Where a MAT MCMC step spends its time (profiling study)

Branch `operator_experiment`, October 2026, M5 Max, one core, Java 21. Purpose: measure the phases of an
MCMC step with the MAT (mutation-annotated tree) likelihood and how they scale with the data, to decide how
(and whether) to multithread it. Full tables: [results.md](results.md). The earlier idea of prefetching
proposals (evaluating K subtree moves in parallel, then accepting the first in sequence) is recorded in the
plan that led here; its ceiling is estimated below.

## Instrumentation (kept, gated, no cost when off)

- `-Dmcmc.profile=<aspects>` selects what `ChainProfile` records (comma-separated; `all` for everything),
  reported after the operator analysis for the chain after the adaptation delay:
  - `steps`: every step split by operator into store, propose, evaluate, accept/restore and adapt (wall
    clock, `MarkovChain`). The operator analysis' `Time(s)`/`Time/Op(us)` columns (whole step, always on)
    come from `operator_study`'s `813361a15`.
  - `likelihoods`: every evaluation of every likelihood component (timed in `CompoundLikelihood` and for
    the posterior itself), with the calculation counts of `Profileable` likelihoods, in total and per
    operator.
  - `mat`: `MATProfile` times the MAT's phases (the resampler's moves, the store's scans, the likelihood's
    evaluation, parsimony, the category Gibbs draw) and counts what they handle, in total and per step of
    the operator in progress; reset with the chain's profile so they cover the same steps.
  The runs below used `-Dmcmc.profile=true -Dmat.profile=true`, the flags of the time; they are now
  `-Dmcmc.profile=steps,mat`.
- Overhead with both on: within noise (-3% to +2% of minutes per million states on four runs).
- `run_profile.sh`, `run_study.sh`, `make_profile_xmls.py` (variants of the benchmark MAT XMLs with the local or
  path block), `simulate_alignment.py` (alignments simulated on big-3's coalescent tree), `summarise.py`
  (the tables), `jfr_top.py` (a flight-recording cross-check).

## Data

| Run | Taxa | Sites | Variable sites | Model | Block |
|---|---|---|---|---|---|
| small-2 | 165 | 1,698 | 294 | JC, strict clock, constant coalescent | local, path |
| wide-1 (VARV) | 64 | 185,578 | 31,181 | JC | local, path |
| big-3 (BDBV) | 899 | 18,903 | 2,405 | HKY | local, path, gamma (4 categories) local |
| big-2 (EBOV) | 1,610 | 18,996 | 18,217 | JC | local, path |
| sim-20k, 100k, 300k | 899 | 20k / 100k / 300k | 2.8k to 135k (14% or 45% of sites) | HKY, simulated on big-3's tree | local |

Chains of 1M states (500k for the path block and the simulations), seed 11, `-Xmx16g`. One seed each;
step times are means over 400k-800k `subtreeLeap` steps, so their noise is far below the differences
discussed. `subtreeJump` in the tables is the `fixedHeightSubtreePruneRegraft` operator (FHSPR).

## Results

### The subtreeLeap step

Microseconds per step. "Resampler" is the MAT hook (history resampling, before + after the move, including
the branch-event updates); "MAT evaluation" is `calculateLikelihood`; "other evaluation" is the rest of the
posterior: the coalescent likelihood and the priors.

| Run | Taxa | Variable sites | STL step | fixed (store+accept/restore+adapt) | STL search+edit | resampler | MAT evaluation | other evaluation (coalescent, priors) | MAT share |
|---|---|---|---|---|---|---|---|---|---|
| small-2-local | 165 | 294 | 17.9 | 1.3 | 1.7 | 7.8 | 2.8 | 4.3 | 59% |
| wide-1-local | 64 | 31181 | 25.9 | 1.0 | 0.9 | 19.6 | 2.5 | 1.8 | 86% |
| big-3-local | 899 | 2405 | 62.0 | 3.5 | 2.8 | 27.0 | 6.7 | 22.0 | 54% |
| big-2-local | 1610 | 18217 | 80.5 | 4.5 | 2.7 | 29.1 | 7.3 | 36.8 | 45% |
| sim-20k-r005 | 899 | 2769 | 42.2 | 3.8 | 1.8 | 9.3 | 6.7 | 20.7 | 38% |
| sim-20k-r02 | 899 | 8988 | 44.7 | 3.7 | 1.5 | 11.1 | 6.8 | 21.6 | 40% |
| sim-100k-r005 | 899 | 13907 | 46.7 | 4.2 | 1.6 | 11.8 | 7.7 | 21.4 | 42% |
| sim-100k-r02 | 899 | 44866 | 53.8 | 4.2 | 1.7 | 17.1 | 8.6 | 22.2 | 48% |
| sim-300k-r02 | 899 | 134581 | 67.7 | 4.5 | 2.0 | 25.9 | 11.8 | 23.6 | 56% |
| big-3-gamma-local | 899 | 2405 | 100.0 | 4.1 | 3.1 | 42.4 | 28.2 | 22.2 | 71% |
| small-2-path | 165 | 294 | 37.2 | 1.7 | 1.8 | 26.2 | 2.9 | 4.5 | 78% |
| wide-1-path | 64 | 31181 | 72.0 | 1.6 | 1.6 | 63.2 | 3.5 | 2.1 | 93% |
| big-3-path | 899 | 2405 | 140.7 | 4.3 | 3.5 | 102.8 | 7.3 | 22.8 | 78% |
| big-2-path | 1610 | 18217 | 171.3 | 5.6 | 4.1 | 116.5 | 8.4 | 36.7 | 73% |

- **The fixed MCMC machinery is small**: 1-6 µs per step (store/restore of all models, adaptation).
- **The coalescent likelihood and priors cost ~22 µs at 899 taxa and ~37 µs at 1,610**, independent of
  the data: about 23 ns per tip, i.e. O(taxa) per tree move (`FastIntervals.calculateIntervals` and
  `CoalescentLikelihood.calculateLogLikelihood`; the flight recordings put them first and second in self time
  on big-3 and big-2). This is a third of a `subtreeLeap` step on big-3 and nearly half on big-2, more than the
  MAT evaluation, and at 10,000 tips it would be ~230 µs per step.
- **The MAT evaluation is 7-12 µs at 899 taxa whatever the number of sites.** It recomputes 14-16 branch
  terms per step (the ~5 changed branches, their children, and the previous rejected step's terms re-dirtied
  on restore) with ~6 transition-matrix cache misses (new branch lengths), plus the O(n) sum (1.7 µs at 899
  taxa, 3.1 at 1,610). With 4 gamma categories it is 28 µs: 45 cache misses per step, each an HKY eigen
  decomposition and 16 logs.
- **The resampler grows slowly with variable sites for subtreeLeap**: 9 µs at 2.8k variable sites to 26 µs
  at 135k (48x the sites, 2.8x the time). The O(sites) scans (`updateBranch` of the reparented branches,
  copy-on-write clones, ambiguous runs) are 1.3 to 6.8 µs of that; most of the growth is in scoring and
  drawing the block's states as events per branch grow. The real data cost more than the simulations at the
  same number of variable sites (big-3: 27 µs at 2.4k) because of ambiguities (hot sites, ambiguity segments).

### Share of the chain's time by component

Over all operators (steps x microseconds per step, chain after the adaptation delay). "MAT likelihood
evaluation" includes the full recomputations forced by the model-parameter moves; "other evaluation" is the
coalescent likelihood and the priors.

| Run | MCMC machinery (store, accept/restore, adapt) | proposal excluding the MAT (tree edit, parameter moves) | MAT history resampling (SPR hooks, Gibbs, categories) | MAT likelihood evaluation | other evaluation (coalescent, priors) | MAT total |
|---|---|---|---|---|---|---|
| small-2-local | 6% | 9% | 37% | 30% | 18% | 67% |
| wide-1-local | 2% | 2% | 89% | 4% | 3% | 94% |
| big-3-local | 5% | 6% | 42% | 19% | 28% | 61% |
| big-2-local | 5% | 5% | 38% | 14% | 38% | 52% |
| sim-20k-r005 | 8% | 6% | 23% | 26% | 38% | 49% |
| sim-100k-r02 | 4% | 3% | 54% | 17% | 21% | 72% |
| sim-300k-r02 | 2% | 2% | 73% | 12% | 11% | 85% |
| big-3-gamma-local | 2% | 2% | 50% | 38% | 9% | 87% |
| big-3-path | 2% | 3% | 76% | 8% | 12% | 83% |
| big-2-path | 2% | 2% | 77% | 5% | 13% | 82% |

### The other operators

Share of the chain's time per operator (steps x microseconds per step):

| Run | subtreeLeap | subtreeJump (FHSPR) | mutationHistoryGibbs | model parameters (rate, kappa, freqs, alpha, upDown) | siteRateCategoryGibbs |
|---|---|---|---|---|---|
| small-2-local | 65% | 8% | 6% | 21% | 0% |
| wide-1-local | 36% | 35% | 28% | 0% | 0% |
| big-3-local | 72% | 12% | 6% | 10% | 0% |
| big-2-local | 75% | 13% | 5% | 6% | 0% |
| sim-20k-r005 | 69% | 12% | 5% | 14% | 0% |
| sim-100k-r02 | 45% | 34% | 12% | 9% | 0% |
| sim-300k-r02 | 28% | 52% | 14% | 6% | 0% |
| big-3-gamma-local | 36% | 5% | 3% | 14% | 42% |
| big-3-path | 64% | 30% | 2% | 4% | 0% |

Scaling with variable sites (899 taxa, local block, microseconds per step of each operator):

| Run | Variable sites | STL step | STL resampler | FHSPR step | FHSPR sites changed per move | Gibbs step | Gibbs sites changed | clock.rate step |
|---|---|---|---|---|---|---|---|---|
| sim-20k-r005 | 2769 | 42.2 | 9.3 | 72.8 | 144 | 28.5 | 14 | 1007 |
| sim-20k-r02 | 8988 | 44.7 | 11.1 | 117.2 | 553 | 49.8 | 55 | 1031 |
| sim-100k-r005 | 13907 | 46.7 | 11.8 | 146.1 | 664 | 61.2 | 69 | 1116 |
| sim-100k-r02 | 44866 | 53.8 | 17.1 | 401.4 | 2642 | 142.9 | 274 | 1290 |
| sim-300k-r02 | 134581 | 67.7 | 25.9 | 1245.9 | 8012 | 346.0 | 834 | 1681 |
| big-3-local | 2405 | 62.0 | 27.0 | 99.1 | 316 | 49.2 | 5 | 1029 |

- **On wide alignments the cost is in FHSPR and the history Gibbs operator, not subtreeLeap.** An FHSPR move
  regrafts at the same height anywhere in the tree, so with the local block the moved subtree changes
  presence along a long chain of ancestors and thousands of sites are redrawn (8,000 per move at 135k variable
  sites; 36 for subtreeLeap): 1.25 ms per step, 52% of the run at 10% of the weight. The Gibbs block update
  redraws every site with an event around its block: 350 µs. Both are O(sites changed) with the sites
  independent given the block.
- **Model-parameter moves recompute every branch term**: 1-1.7 ms per step (1,797 terms, one cache miss
  each), 3 ms with gamma; 6-14% of the run at their usual weights.
- **With gamma categories the category Gibbs operator is 42% of the run**: 11 ms per step (8.3 ms for the
  draw over all sites and nodes, 3 ms for the evaluation it forces), at weight 8.99.
- The path block costs 3.5-4x the local block for subtreeLeap and 6x for FHSPR, as before.

## What this says about multithreading

Ceilings for a `subtreeLeap`-dominated run on big-3 or sim-300k, with 4 threads, from the tables (Amdahl,
no thread hand-off cost, which at these step lengths would take a further 10-20%):

| Approach | What it parallelises | Serial residue | Ceiling |
|---|---|---|---|
| A. Prefetching proposals | the resampler + MAT evaluation of each proposal (34-38 µs of 62-68) | fixed + search + coalescent/priors + applying the delta: 28-30 µs per assessed proposal; 2.8 of 4 proposals useful | ~1.5x |
| B. Site-range splitting inside a move | for subtreeLeap only the O(sites) scans (2-7 µs): nothing to gain; for FHSPR and Gibbs the per-site draws (900 of 1,250 µs at 135k sites) | the rest of those moves and all other steps | 1.4x on sim-300k (FHSPR is 52% of the run), nothing on big-3 |
| C. Branch-parallel evaluation | the model-parameter moves (1-3 ms each, 6-14%), and the category Gibbs draw (42% of a gamma run) | everything else | ~1.1x; ~1.6x on gamma runs |
| E. Independent chains | everything | burn-in | linear |

No threading scheme reaches 2x on these data, because the step is a sum of several comparable parts and the
largest non-MAT part, the coalescent, is serial and grows with the number of tips. The serial fixes the
profile points at are larger and apply to every operator:

1. **An incremental coalescent likelihood.** Recomputing all intervals and the density after every tree move
   is 22 µs at 899 tips and 37 at 1,610 (35-45% of a subtreeLeap step), and O(tips). A subtreeLeap moves one
   node height: the interval list changes locally and the constant-size density is a sum over intervals, so
   an O(changed intervals) update is possible (and `FastIntervals` already sorts incrementally).
2. **Cheaper branch terms.** Half the MAT evaluation is transition-matrix cache misses (an HKY eigen
   decomposition plus 16 logs per new branch length; 45 per step with 4 categories): closed-form HKY/JC
   probabilities, or per-branch matrices keyed on the branch rather than the length, would cut the gamma
   evaluation two- to three-fold. Fewer terms per step (the re-dirtying of the previous rejected step's
   terms on restore roughly doubles them) would help all runs.
3. **FHSPR and the Gibbs block on wide alignments.** Their cost is algorithmic: thousands of sites redrawn per
   move. Options are a lower weight relative to subtreeLeap, a block that does not grow with the regraft
   distance, or, as the only place where threads clearly fit, drawing the independent sites of a block in
   parallel (B), which needs per-thread random streams and per-thread writes to the store.
4. **The category Gibbs draw** (8 ms over all sites and nodes) is the other site-parallel candidate, and the
   first thing to attack in gamma runs.

Independent chains remain the best use of several cores for the MAT.

## Caveats

- The estimate at the start of this study, that the fixed MCMC machinery was ~50 µs per step, came from the
  operator-study runs made with an older frozen jar; on this branch it is 1-6 µs. The timings here are from
  the current build.
- The flight recordings agree on the ranking (coalescent first, then the resampler's block sampling, then the
  MAT branch terms) but attribute a larger share to the coalescent (46-64% of all samples) than the timers
  do (30-45%); JFR collected far fewer samples than expected on this machine, so the timers are taken as
  authoritative.
- The simulated alignments have no ambiguities; the real alignments' ambiguity runs make the local block
  about twice as expensive at the same number of variable sites (big-3 against sim-20k).
- Single seed per run; the chain starts at the simulation's tree for the simulated data and at a coalescent
  tree for the real data, with 1-5% of the steps as adaptation delay excluded from the step profile.

## Reproducing

```
ant dist
python3 make_profile_xmls.py ../../../beast-benchmark/Benchmarks xmls
# the big-3 coalescent tree: a <report><treeModel idref="treeModel"/></report> on the taxa + coalescentSimulator of the XML
python3 simulate_alignment.py ../../../beast-benchmark/Benchmarks/Benchmark-big/benchmark-big-3-mat.xml big3.tree xmls/sim-100k-r02.xml --sites 100000 --rate 0.02 --chain-length 500000
./run_study.sh xmls runs          # profiled, unprofiled and JFR groups; JAVA_OPTS=-Xmx16g
python3 summarise.py runs --markdown results.md
JAVA_HOME=$(/usr/libexec/java_home -v 21) python3 jfr_top.py runs/jfr/big-3-local/run.jfr
```
