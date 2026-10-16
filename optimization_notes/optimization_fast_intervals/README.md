# optimization_fast_intervals

Branch: `optimization_fast_intervals` (from `main`)

| Commit | Change |
|---|---|
| `a904f313e` | Sort coalescent interval times incrementally in FastIntervals |
| `d0e9ad8bd` | Copy event times and finish time in `FastIntervals.copyIntervals` (bug fix) |

Only `src/dr/evolution/coalescent/FastIntervals.java` changes, plus its test
`src/test/dr/evolution/coalescent/FastIntervalsTest.java`. `TreeIntervals`, the node-mapping
`Intervals` class and the code that uses node mapping are untouched. FastIntervals is only used for
straightforward coalescent models: `<treeIntervals>` elements, which BEAUti writes for the parametric
coalescents (constant, exponential, logistic, expansion), and `<coalescentLikelihood>` with
`<populationTree>`.

## Why

Profiling MCMC runs outside the data likelihood showed that calculating the coalescent intervals took
11–16% of the run time on large trees. Most of that (9–12%) was the two `Arrays.sort` calls in
`FastIntervals.calculateIntervals()`, which fully sort the sample and coalescent times after every
change to the tree. That is O(n log n), even though most moves change only one or a few node heights.

## 1. Incremental sorting (`a904f313e`)

`TreeIntervals` adds events by traversing the nodes in node-number order, using the
`addSampleEvent(time, nodeNumber)` and `addCoalescentEvent(time, nodeNumber)` methods. So the time at
array index *i* always belongs to the same node, and between evaluations only the entries for nodes
whose heights changed differ.

FastIntervals now keeps, for each of the sample and coalescent times, an `IncrementalSorter` that
remembers the previous times and their sorted order (as indices). On the next calculation it:

1. finds the indices whose times changed (*k* of them), by comparing with the previous times;
2. sorts those *k* indices by time (insertion sort for small *k*, otherwise a stable merge sort);
3. merges them with the unchanged indices, which are already in order, in one pass;
4. gathers the sorted times.

This costs O(n + k log k) rather than O(n log n). Ties are ordered exactly as a stable sort would
order them, so the sorted times are identical to `Arrays.sort`, and so are the intervals.

Details:

- If events are added without node numbers (`addSampleEvent(time)` / `addCoalescentEvent(time)`), the
  order is not guaranteed, so FastIntervals falls back to `Arrays.sort` and invalidates the cached
  order.
- `copyIntervals` also copies the sorter state, so after a restore (`TreeIntervals` swaps the stored
  copy back in) the next calculation is still incremental.
- When many times change (e.g. scaling all node heights), the merge still works. It just costs close
  to a full sort.

## 2. Bug fix: `copyIntervals` (`d0e9ad8bd`)

This bug is already on `main`. `copyIntervals` copied the interval lengths, types and lineage counts
and the start time, but not `eventTimes` or `finishTime`. After a rejected move, `TreeIntervals`
swaps the stored copy back in. That copy then returned stale values from `getIntervalTime(i)` and
`getTotalDuration()`, left over from whenever that instance last calculated its own intervals
(possibly thousands of states earlier). The values stayed stale until the next move that changed the
tree. Also, if the source's intervals had not yet been calculated, the copy recalculated from its own
stale events. It now copies the source's events in that case.

Consequences for real runs are very limited:

- `CoalescentLikelihood` (BEAUti's parametric coalescents) reads only interval lengths, types, lineage
  counts and the start time. These were all copied correctly, so it was not affected.
- BEAUti skygrid and skyride analyses use `<populationTree>`, which builds node-mapped `Intervals`
  rather than `FastIntervals`, so they were not affected. SkyGlide, SmoothSkygrid and the time-varying
  clock use `BigFastTreeIntervals`, so they were not affected either.
- Hand-written XMLs that give a `<treeIntervals>` (FastIntervals) to a skygrid or skyride via
  `<intervals>` are affected in one case: with HMC on node heights using the skygrid gradient. The
  likelihoods themselves only read times after a tree change, which recalculates first, so they are
  fine. But `GMRFGradient` reads `getIntervalTime` at the current state, before the tree moves.
  Straight after a rejection it would get the wrong gradient at the start of the trajectory, which
  breaks HMC's reversibility (a possible small bias, as well as worse mixing). Any statistic that logs
  interval times could also have logged stale values.

The fix costs one extra copy of `eventTimes` per `storeState`.

## Testing

`FastIntervalsTest` (JUnit):

- `testIncrementalSortingMatchesFullSort`: 5,000 random updates of none, one, a few or many node
  heights (including ties and tip-time changes), with store/restore cycles done as `TreeIntervals`
  does them. After each update, every interval length, time, type and lineage count must equal a full
  sort exactly. The restored copy is also checked before its events are added again; this check fails
  without the `copyIntervals` fix.
- `testEventsWithoutNodeNumbers`: events added without node numbers, in any order, are fully sorted.
- `testCopyBeforeCalculation`: a copy taken before the intervals are calculated matches its source;
  this also fails without the fix.

MCMC logs with the same seed are bit-identical (apart from the header lines) between `main` and this
branch, for the test XML here and for the big-2 and big-3 benchmark baselines.

## Performance

Machine: Apple M5 Max, OpenJDK 26.0.1, single-threaded, one run at a time.

### Coalescent-only benchmark (`prior-only-1610.xml`)

This file is the 1,610 dated EBOV tips from `benchmark-big-2`, with the tree sampled from a
constant-size coalescent alone: no alignment, substitution model, clock or tree likelihood. That makes
the coalescent and tree machinery most of the cost. Operators: subtreeLeap (weight 1610),
fixedHeightSPR (161), uniform on internal node heights (161), up on node heights and popSize (3), and
scale on popSize (3). popSize has a 1/x prior. 10,000,000 states, seed 1234.

| FastIntervals | Run time | min / million states | Speed-up |
|---|---|---|---|
| `main` | 19.10 min | 1.91 | – |
| incremental sort (`a904f313e`) | 13.07 min | 1.31 | **31.5% less time** (1.46×) |
| + `copyIntervals` fix (`d0e9ad8bd`) | (running) | | |

The final log lines are identical (joint −677.9039, root height 2.15861 at state 10,000,000).

To run:

```
beast -seed 1234 -overwrite prior-only-1610.xml
```

### Microbenchmark (cost per evaluation, 1,610 tips, one node height changed)

| | `main` | incremental |
|---|---|---|
| `calculateIntervals` (sorting and building the intervals) | ≈ 51 µs | ≈ 12 µs |

The remaining cost per evaluation with the incremental sort breaks down as follows:

| Step | Cost |
|---|---|
| `copyIntervals` (store) | 1.25 µs |
| adding events (done by `TreeIntervals`, not changed) | 2.74 µs |
| `calculateIntervals` (change detection, merge, interval loop) | 7.34 µs |

### With a data likelihood

On `benchmark-big-2` (1,610 taxa) with the single-threaded MAT likelihood, the run was 6.5% faster.
The likelihood dominates there, so the saving is a smaller share of the total. BEAGLE multi-threaded
timings vary by about ±10% from run to run, which is too noisy to measure this change.
