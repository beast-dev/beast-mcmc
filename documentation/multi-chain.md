# Running multiple chains in BEAST

BEAST can run several MCMC chains at once in one process, on separate threads:

- **Independent chains** (`-chains K`): K chains from different random starting points,
  each writing its own output files. While they run, BEAST compares them to estimate a
  common burn-in and to report how well they agree.
- **Parallel tempering**, also called Metropolis-coupled MCMC or MC³ (`-mc3_chains N`): a
  set of one cold chain and N−1 heated chains that occasionally swap states. This helps
  the cold chain escape local optima. Only the cold chain's samples are logged.

The two combine. `-chains 4 -mc3_chains 3` runs four independent parallel-tempering sets,
12 chains in total.

With neither option, BEAST runs a single chain exactly as it always has: same output, same
file names.

---

## Quick start

```bash
# four independent chains
beast -chains 4 run.xml

# one parallel tempering set: a cold chain and three heated chains
beast -mc3_chains 4 run.xml

# two independent tempering sets of three chains each, at most 4 chains computing at once
beast -chains 2 -mc3_chains 3 -chain_threads 4 run.xml
```

### Options

| Option | Default | Meaning |
|---|---|---|
| `-chains K` | 1 | Number of independent chains (or independent tempering sets). |
| `-chain_threads T` | one per processor | Maximum number of chains computing at the same time. |
| `-convergence_every S` | 60 | Seconds between convergence reports when K > 1 (0 for none). |
| `-mc3_chains N` | 1 | Chains in each tempering set: one cold chain plus N−1 heated ones. |
| `-mc3_delta d` | 1.0 | Heating increment: chain *i* gets temperature 1/(1 + d·*i*). |
| `-mc3_temperatures t1,t2,…` | | The heated chains' temperatures, instead of `-mc3_delta` (N−1 values). |
| `-mc3_swap n` | 100 | Attempt temperature swaps every *n* states. |
| `-mc3_scheme name` | `original` | Swap scheme: `original`, `stochastic_single`, `stochastic_multiple`, `deterministic_single`, `deterministic_multiple` (see below). |
| `-seed s` | clock | Seed for the whole run. Each chain's own seed is derived from it. |
| `-thread_scheme s` | `auto` with several chains | How each chain uses its threads: `serial`, `split` or `auto` (see [Threading](#threading)). |
| `-threads n` | set by the scheme | Threads for each chain's likelihood calculations (overrides the scheme). |
| `-beagle_threads n` | set by the scheme | CPU threads for each BEAGLE instance (overrides the scheme). |

---

## Independent chains

Each chain is built from the same XML file, but each has:

- **its own random number stream.** The seeds are drawn from `-seed` and printed at the
  start, so chains start from different random states (for example, different random
  starting trees).
- **its own output files.** A `.chainK` suffix goes before the extension:
  `run.chain1.log`, `run.chain1.trees`, `run.chain2.log`, and so on. Checkpoint files get
  the same suffix.
- **its own model objects.** Nothing is shared between chains.

Only chain 1 prints its screen log; the other chains report through the convergence
monitor. At the end, each chain's operator analysis is printed under a `Chain K:` heading.

Elements after `<mcmc>` in the XML, such as `<report>` or a marginal likelihood estimator,
are not run in multi-chain mode.

### Convergence monitor

When K > 1, BEAST keeps each chain's samples in memory. It compares the chains every
`-convergence_every` seconds and once more at the end:

```
Convergence at state 598000 (3 chains):
  Chain 1    state 612000       posterior -192975.1234
  Chain 2    state 604000       posterior -192981.5678
  Chain 3    state 598000       posterior -192969.9012
  common burn-in 108000 states (18%): posterior R-hat 1.0094
  min ESS 412.3 (frequencies1)   max R-hat 1.0151 (ag)   ASDSF 0.0031 (64 clades)
```

**Common burn-in.** All chains share one burn-in, chosen by when they agree rather than by
looking at each trace separately. Burn-ins from 0% to 50% of the samples are tried in 2%
steps. The burn-in chosen is the smallest one after which the split R-hat of the chains'
posteriors stays below 1.01, for that burn-in and for every larger one tried.

The "every larger one" condition matters. The initial climb inflates each chain's own
variance, so R-hat with no burn-in can look fine while the chains still disagree. If no
burn-in up to 50% qualifies, the report says so:

```
  no common burn-in yet: posterior R-hat 1.1613 after 50% (needs < 1.01); below is after 50%
```

That means the chains need to run longer, not that more should be discarded.

**After the burn-in**, the monitor reports three diagnostics:

| Diagnostic | Meaning | Rule of thumb |
|---|---|---|
| min ESS | Smallest effective sample size of any logged value, summed over the chains. | > 200 |
| max R-hat | Largest split R-hat (Gelman et al., *BDA3*) of any logged value. Each chain is split in half, so a trend within a chain shows up as well as disagreement between chains. | < 1.01 |
| ASDSF | Average standard deviation of clade frequencies across the chains, over clades with frequency ≥ 0.1 in some chain. | < 0.01 |

The values used are those in each chain's main log file (the logger with the most numeric
columns), at its logging interval. The clades come from the first tree logger. Trees are
not kept in memory: each clade is reduced to a 64-bit hash and counted. This keeps memory
use small even for trees with thousands of tips.

The monitor only reports; it never stops the chains. The same figures are written, one row
per report, to `<xml name>.convergence.log`:

```
state   burnin  posteriorRhat  minESS  minESS.parameter  maxRhat  maxRhat.parameter  ASDSF   clades
598000  108000  1.0094         412.3   frequencies1      1.0151   ag                 0.0031  64
```

### Checkpointing

`-save_every`, `-save_at`, `-save_state` and `-load_state` work per chain. With
`-save_state run.state`, chain K saves to `run.chainK.state`. With `-load_state run.state`,
chain K loads `run.chainK.state`. Each checkpoint contains only its own chain's state,
including its random number stream.

Checkpointing is not yet available with parallel tempering (`-mc3_chains` > 1). BEAST stops
with an error if you ask for both.

### Reproducibility

The same `-seed` gives the same output for each chain, whatever the number of chain
threads. One caveat: changing the number of threads used inside a chain's likelihood
calculation (`-threads`, `-beagle_threads`) changes the order of floating-point sums, and
chains then drift apart in the last digits. By default these thread counts are set by the
threading scheme from each chain's thread budget, which depends on `-chain_threads` and the
number of processors (see [Threading](#threading)). For bit-for-bit repeatable runs on
different machines or with different `-chain_threads`, give `-threads` and `-beagle_threads`
explicitly.

---

## Parallel tempering (MC³)

A tempering set has one cold chain, which samples the posterior, and N−1 heated chains.
These sample flattened versions of it, so they cross valleys in the posterior more easily.

Each chain's "temperature" in BEAST is the power applied to the posterior: 1.0 for the
cold chain, smaller for heated chains. With `-mc3_delta d`, chain *i* uses 1/(1 + d·*i*):
with d = 1.0 that is 1, 0.5, 0.33, 0.25, …. A smaller delta gives gentler heating and more
accepted swaps.

All chains in a set run in parallel for `-mc3_swap` states. They then wait for each other
while BEAST proposes swapping the temperatures of pairs of chains, accepting with the usual
Metropolis probability. When a swap involves the cold chain, the cold chain's log files
move with the temperature. The output therefore always comes from whichever chain is cold
at the time, and only one set of files is written. The operator tunings move with the
temperatures too.

Swap schemes (`-mc3_scheme`):

| Scheme | Which pairs are proposed at each swap |
|---|---|
| `original` | One pair of chains chosen at random (any two temperatures). |
| `deterministic_single` / `deterministic_multiple` | Neighbouring temperatures only, alternating between the pairs (0,1),(2,3),… and (1,2),(3,4),… One random pair, or all of them. |
| `stochastic_single` / `stochastic_multiple` | As above, but the alternation is chosen at random each time. |

At the end, a table of attempted and accepted swaps for each pair of temperatures is
printed. Acceptance between neighbouring temperatures of roughly 20–60% is a reasonable
target. If it is near zero, reduce `-mc3_delta`. If it is near one, the heated chains are
too cold to help.

Each chain in a set costs as much to compute as the cold chain. A set of N chains needs N
times the processing of a single chain, and the chains wait for the slowest one at each
swap.

---

## Threading

### How the work is divided

There are three levels of parallelism, from outermost to innermost:

1. **Chains.** Every chain, cold or heated, has its own dedicated thread. A chain is always
   computed on the same thread, which BEAGLE GPU instances need. At most `-chain_threads`
   chains compute at any moment; the rest wait their turn. Independent chains run in
   segments of 1% of the chain length, so with fewer chain threads than chains they take
   turns and progress together. Members of a tempering set synchronise at every swap
   interval. Independent sets never wait for each other.
2. **Likelihood components** (`-threads`). A chain's likelihoods (one per data partition)
   and priors can be computed in parallel by a thread pool within the chain.
3. **BEAGLE CPU threads** (`-beagle_threads`). Each BEAGLE CPU instance splits its site
   patterns over threads. There is one instance per data partition (or one for several
   partitions with `-beagle_multipartition`). **This count is per instance, not per chain.**

Levels 2 and 3 multiply: a chain computing 3 partitions in parallel, each with 4 BEAGLE
threads, keeps 12 threads busy. Busier than the processors allow, and threads that must
finish together wait for each other, so the inner levels should not both be given the whole
budget.

### Threading schemes

Each chain has a **thread budget**:
- **One chain:** the number of processors.
- **Several chains:** the processors shared out between the chains computing at once,
  max(1, P / `-chain_threads`), where P is the number of processors.

The scheme (`-thread_scheme`) decides how the chain spends its budget on levels 2 and 3:

| Scheme | Likelihood components | BEAGLE CPU threads per instance |
|---|---|---|
| `serial` | one after another | the whole budget |
| `split` | in parallel (one thread per BEAGLE instance) | the budget divided between the instances |
| `auto` | `split` or `serial`, chosen for each chain from its data | as for the scheme chosen |

`auto` chooses `split` when either of these holds; otherwise it chooses `serial`:
- the chain has at least two BEAGLE GPU instances;
- the chain has at least two BEAGLE CPU instances that each have enough work, and a budget
  of at least two threads.

"Enough work" means at least 100,000, measured as tips × patterns × (states / 4)², so
1000 tips and 100 nucleotide patterns is just enough. The choice is made when each chain
starts, so it is based on the actual BEAGLE instances and where they run. The chosen scheme
is printed, for example:

```
Threading scheme: split (auto), 9 threads per chain
  likelihoods computed in parallel (3 threads)
  BEAGLE: 3 CPU instances with 3 threads each, 0 GPU instances
```

**Defaults:**
- **Several chains:** `auto`.
- **One chain:** BEAST's usual threading (one thread per likelihood component, and one BEAGLE
  thread per processor for every instance) unless `-thread_scheme` is given.

`-threads` or `-beagle_threads`, if given, override that level and the scheme sets only the
other.

P is the number of processors Java reports, which counts hyper-threads. On a machine with
hyper-threading, P is usually twice the number of physical cores. On a machine whose cores
differ in speed (for example Apple's performance and efficiency cores), it counts all of them;
see [Benchmarks](#benchmarks).

### Memory

Each chain holds a complete copy of the model, including the BEAGLE partial-likelihood
buffers. A run with K independent sets of N tempered chains needs roughly K × N times the
memory of a single chain. Increase the Java heap (`-Xmx`) accordingly. GPU memory is used in
the same proportion if all chains are on one GPU.

---

## Tuning for your hardware

The aim is to keep every core (or GPU) busy without oversubscribing it. Rules of thumb:

- **Total CPU threads.** Chains computing at once × threads per chain should be about the
  number of physical cores. The `serial` and `split` schemes keep each chain within its
  budget. Setting `-threads` and `-beagle_threads` yourself can exceed it.
- **Chains before threads.** For the same cores, several chains with a few threads each
  compute more states in total than one chain with all the threads (2–2.4× in the
  benchmarks). But each chain is slower and has its own burn-in, so choose the number of
  chains you want for the analysis first, then give them the cores.
- **Small data, few patterns.** Threads within a chain gain little and computing partitions
  in parallel can be much slower. Use more chains with 1–2 threads each.
- **Several large partitions.** Computing them in parallel (`split`) helps most, especially
  with few chains.
- **Hyper-threading** rarely helps the likelihood calculations. Count physical cores, and on
  machines with performance and efficiency cores, count the performance cores.

### CPU only

**Laptop, 8 physical cores (16 hyper-threads), 4 independent chains:**

```bash
beast -chains 4 -chain_threads 4 -beagle_threads 2 run.xml
```

All four chains run at once with 2 BEAGLE threads per instance. Left to the defaults, each
chain would have a budget of 16 / 4 = 4 threads, counting hyper-threads: with one partition
(`serial`), all 4 would go to BEAGLE, and with several sizeable partitions (`split`), they
would be divided between them. Giving `-beagle_threads` yourself fixes the BEAGLE threads
per instance; with several partitions computed in parallel that is 2 threads for each.

**Workstation, 32 cores, 4 tempering sets of 4 chains (16 chains):**

```bash
beast -chains 4 -mc3_chains 4 -chain_threads 16 -beagle_threads 2 run.xml
```

**More chains than cores**, for example 8 chains on 4 cores, for a quick look at many
independent starts:

```bash
beast -chains 8 -chain_threads 4 -beagle_threads 1 run.xml
```

The chains take turns in segments and all advance together. Total throughput is the same as
running them in two batches, but the convergence report compares all eight throughout.

**Several partitions:** each partition is a separate BEAGLE instance. With the `split`
scheme (which `auto` chooses for several sizeable partitions), the chain's budget is divided
between them and they are computed in parallel. If you set `-beagle_threads` yourself, each
instance gets that many threads: with 3 partitions and `-beagle_threads 4`, a chain computing
them in parallel uses 12 threads.

### GPUs

Use `beast -beagle_info` to list the BEAGLE resources and their numbers. Resource 0 is the
CPU; GPUs are usually 1, 2, ….

BEAGLE instances are created in a fixed order: chain 1 first (its cold chain, then each
heated chain), then chain 2, and so on, with each chain's partitions in XML order.
`-beagle_order` gives the resource for each instance in turn. When the list runs out it
starts again from the beginning. This is how chains are spread over GPUs.

**One GPU, 2 independent chains, one partition:**

```bash
beast -beagle_GPU -chains 2 run.xml
```

Both chains share the GPU. This works well when one chain cannot keep the GPU busy, which
is common for moderate data sets.

**Two GPUs, 2 independent chains, one partition each:**

```bash
beast -beagle_GPU -beagle_order 1,2 -chains 2 run.xml
```

Chain 1's instance goes to GPU 1 and chain 2's to GPU 2.

**Two GPUs, 2 chains, 2 partitions each** (4 instances):

```bash
beast -beagle_GPU -beagle_order 1,1,2,2 -chains 2 run.xml
```

**One GPU and the CPU, 2 chains, one partition:**

```bash
beast -beagle_order 1,0 -chains 2 run.xml
```

Chain 1 goes on the GPU and chain 2 on the CPU. The chains will run at different speeds; the
convergence report compares them at the state the slower one has reached.

**Parallel tempering on GPUs:** each heated chain is a separate instance. With one partition,
`-mc3_chains 4 -beagle_order 1,1,2,2` puts the cold chain and the first heated chain on GPU 1,
and the other two heated chains on GPU 2.

A chain running on a GPU uses little CPU; its thread mostly waits for the GPU. You can
therefore allow as many chain threads as there are chains (the default when there are
enough processors), even with fewer cores than chains.

### Checking the result

At the start, BEAST prints which BEAGLE resource each of chain 1's instances is using
(`Using BEAGLE resource 1: …`). Then, as the first chain starts, it prints the threading
scheme it chose, which is what the chains actually use (the BEAGLE thread counts printed
while the model is built are replaced by the scheme's). The setup messages of the other
chains are not shown, because they would repeat chain 1's. The other chains follow the same
rules, so their resources can be worked out from `-beagle_order` as described above.

The screen log reports the speed of chain 1 as time per million states. Compare it with a
single-chain run to see how much running chains side by side slows each one down.

---

## Benchmarks

The threading advice above comes from benchmarks on one machine: an Apple M5 Max with 18
cores (6 performance and 12 efficiency cores), BEAGLE 4 on the CPU, October 2026.

**Datasets.** Most results use `benchmark-big-3` from the BEAST benchmark set (899 taxa,
18,900 sites, HKY, coalescent), split into:
- **p1:** 1 partition;
- **p3:** 3 partitions, one per codon position;
- **p12:** 12 partitions, blocks of sites.

Smaller datasets were also used:
- **small:** `benchmark-small-2`, 165 taxa and 258 patterns;
- **tall:** `benchmark-tall-1`, 1441 taxa and 593 patterns.

**Measure.** Throughput is states per second, summed over the chains, with start-up time
subtracted. Each figure is the mean of 2–3 runs; runs varied by up to about 15%.

**Computing partitions in parallel.** Throughput relative to computing them one after another
(`serial`), for the same thread budget per chain. With 1 chain the budget was the 6
performance cores; with 4 chains it was 4 threads each.

| Dataset | 1 chain: split | 1 chain: one thread per partition | 1 chain: both levels at full budget | 4 chains: split | 4 chains: both at full budget |
|---|---|---|---|---|---|
| big-3, 1 partition | – | – | 0.94 | – | 0.93 |
| big-3, 3 partitions | **1.15** | 0.93 | 1.02 | 1.08 | 1.08 |
| big-3, 12 partitions | **1.75** | 1.69 | 1.54 | **1.52** | 1.23 |
| tall, 3 partitions | **1.22** | 1.12 | 1.19 | 0.99 | 0.86 |
| small, 3 partitions | 0.68 | 0.69 | 0.69 | 0.51 | 0.50 |

What this shows:
- **Parallel partitions pay off only with enough data per partition.** They gained up to 75%
  for large partitions, cost 30–50% for small ones, and were slightly slower with one
  partition.
- **Split is the best way to compute partitions in parallel.**
- **Giving both levels the whole budget (the old multi-chain default) was never best.**

**All 18 cores for one chain** was 14–29% slower than the 6 performance cores. BEAGLE's
threads wait for the slower efficiency cores.

**Chains versus threads.** Total throughput relative to one chain using all the threads:

| Cores | big-3, 1 partition | big-3, 3 partitions |
|---|---|---|
| 18 | 18 chains × 1 thread: 2.1× (3 chains × 6: 1.6×) | 3 chains: 2.4× |
| 6 performance cores | 2–3 chains: 1.3× | 3 chains: 2.0× |

**The new defaults against the old.** States per second, relative to the old behaviour on
the same runs:
- **4 chains:** `auto` against the old default of `-threads` and `-beagle_threads` both
  set to the per-chain budget.
- **1 chain:** `-thread_scheme auto` against BEAST's usual threading.

| Dataset | 4 chains | 1 chain |
|---|---|---|
| big-3, 1 partition | 1.07 | 0.93 |
| big-3, 3 partitions | 1.11 | 1.27 |
| big-3, 12 partitions | 1.19 | 1.03 |
| small, 3 partitions | 1.84 | 1.35 |

**Caveats.** Machines with uniform cores, more memory bandwidth or different BEAGLE builds
may differ, so it is worth timing a short run of your own analysis. The GPU parts of the
`auto` rule could not be tested on this machine.
