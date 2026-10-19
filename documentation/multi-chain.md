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
| `-threads n` | see [Threading](#threading) | Threads for each chain's likelihood calculations. |
| `-beagle_threads n` | see [Threading](#threading) | CPU threads for each BEAGLE instance. |

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
chains then drift apart in the last digits. By default these thread counts depend on
`-chain_threads` and on the number of processors (see below). For bit-for-bit repeatable
runs on different machines or with different `-chain_threads`, give `-threads` and
`-beagle_threads` explicitly.

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
2. **Likelihood components** (`-threads`). A model's likelihoods and priors can be
   evaluated in parallel by a thread pool within each chain.
3. **BEAGLE CPU threads** (`-beagle_threads`). Each BEAGLE CPU instance splits its site
   patterns over threads. There is one instance per data partition (or one for several
   partitions with `-beagle_multipartition`). **This count is per instance, not per chain.**

### Defaults

The defaults depend on whether one chain or several are running.

**One chain:**
- `-threads`: one thread per likelihood component.
- `-beagle_threads`: one per processor, for every BEAGLE instance.

**Several chains:** the processors are shared between the chains computing at once. With
P processors and T = `-chain_threads`:
- `-chain_threads`: min(P, total number of chains).
- `-threads`: max(1, P / T). When this is 1, components are computed one after another.
- `-beagle_threads`: max(1, P / T).

Any option you give explicitly is used as given.

P is the number of processors Java reports, which counts hyper-threads. On a machine with
hyper-threading, P is usually twice the number of physical cores.

### Memory

Each chain holds a complete copy of the model, including the BEAGLE partial-likelihood
buffers. A run with K independent sets of N tempered chains needs roughly K × N times the
memory of a single chain. Increase the Java heap (`-Xmx`) accordingly. GPU memory is used in
the same proportion if all chains are on one GPU.

---

## Tuning for your hardware

The aim is to keep every core (or GPU) busy without oversubscribing it. Rules of thumb:

- **Total CPU threads.** Chains computing at once × threads per chain should be about the
  number of physical cores. Threads per chain is roughly `-beagle_threads` × the number of
  BEAGLE instances computing at the same time.
- **Small data, few patterns.** Chains scale much better than threads within a chain.
  Prefer more chains with 1–2 threads each.
- **Large data, many patterns.** BEAGLE threading scales well, so fewer chains with more
  threads each is fine.
- **Hyper-threading** rarely helps the likelihood calculations. Count physical cores.

### CPU only

**Laptop, 8 physical cores (16 hyper-threads), 4 independent chains:**

```bash
beast -chains 4 -chain_threads 4 -beagle_threads 2 run.xml
```

All four chains run at once with 2 BEAGLE threads each. Left to the defaults, each chain
would get 16 / 4 = 4 threads, counting hyper-threads.

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

**Several partitions:** each partition is a separate BEAGLE instance, and each gets
`-beagle_threads` threads. With 3 partitions and `-beagle_threads 4`, a chain can use 12
threads. Reduce `-beagle_threads`, or use `-threads 1` so the partitions are computed one
after another, to stay within the core count.

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

At the start, BEAST prints which BEAGLE resource and how many CPU threads each of chain 1's
instances is using. Look for lines like `Using BEAGLE resource 1: …` and `Using 2 threads
for CPU.`. The setup messages of the other chains are not shown, because they would repeat
chain 1's. The other chains follow the same rules, so their resources can be worked out from
`-beagle_order` as described above.

The screen log reports the speed of chain 1 as time per million states. Compare it with a
single-chain run to see how much running chains side by side slows each one down.
