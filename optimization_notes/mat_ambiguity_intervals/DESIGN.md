# MAT: ambiguities as intervals (design)

Branch: `operator_experiment`. Status: design, for agreement before implementation.

## Why

Delphy represents missing data on each branch as intervals of sites ("missations"), not one entry per
site. The BEAST MAT stores one ambiguity event per site per branch. On the BDBV benchmark (big-3: 899
genomes, 2,405 variable sites, many runs of N), at initialisation:

| Branches | Ambiguity events | Runs of consecutive sites | Events per branch (mean / max) | Substitutions |
|---|---|---|---|---|
| 899 tip | 283,785 | 2,517 | 316 / 1,512 | 3,909 |
| 897 internal | 71,794 | 396 | 80 / 1,379 | 331 |

So there are 355,579 ambiguity events against a few thousand substitutions, but only 2,913 runs (122
events per run). Everything that loops over a branch's events, or over the sites that have ambiguity
events around a resampling block, pays for every N individually.

Profile of big-3 integrated hits (`operator_experiment` at `30911acba`, fast tree model, 300,000
states):

| | Share of run time | Loops over ambiguities per site? |
|---|---|---|
| History resampler | 70% | |
| &nbsp;&nbsp;building per-group conditionals (`Conditional`, `UniformConditional`) | 36% | no (see "Related" below) |
| &nbsp;&nbsp;per-site handling of ambiguity-only sites (`sampleAmbiguitySites` itself) | 12% | yes: about 1,150 sites per move |
| &nbsp;&nbsp;finding candidate sites, writing ambiguity flags, updating branch events | about 10% | yes |
| MAT likelihood | 17% | |
| &nbsp;&nbsp;`getComposition` (and `getExitRate`) | 9% | yes: about 13 nodes per evaluation, 46 events each |
| &nbsp;&nbsp;branch log likelihoods | 4% | partly |
| Coalescent, FastIntervals, FastTreeModel, store/restore | about 12% | no |

About 30% of the run time is spent on per-site ambiguity work that intervals would reduce to per-run
work.

## Design

### 1. Representation of ambiguity events

`BranchMutations` keeps substitutions as now (a list sorted by pattern, with optional hit times), and
replaces the per-site ambiguity entries by intervals:

- `int[] ambiguityStarts, ambiguityEnds`: half-open pattern ranges `[start, end)`, sorted and disjoint,
  of the patterns at which the node is ambiguous and its parent is not;
- `getAmbiguityIntervalCount()`, `getAmbiguityStart(i)`, `getAmbiguityEnd(i)`,
  `getAmbiguityCount()` (the number of patterns), and `isAmbiguous(p)` (binary search);
- the per-site accessors `isAmbiguity(k)` and the ambiguity entries in `size()`/`getPattern(k)` go;
  callers that loop over all events are changed to loop over substitutions and intervals separately.

The state of the parent at an ambiguous pattern (the "from" state) is not stored: it is the parent's
state, which the store has.

Runs only compress if consecutive patterns are consecutive sites. That holds for the compressions the
MAT uses (uncompressed, `constantOnly`, `ambiguousConstant`: variable patterns are kept in site order).
With other compressions intervals are still correct, just not shorter.

### 2. Node states stay dense

The store keeps one byte per pattern per node (`-1` for ambiguous) as the source of truth. The
resampler and parsimony read and write them, and they make the change incremental: intervals become a
derived view of the events on each branch, like the substitution lists are now. Replacing the dense
states by a reference sequence plus per-node differences (as Delphy does) would remove the remaining
O(sites) scans but touches everything; it is out of scope.

### 3. Building the events

- `updateBranch` (a branch whose parent changed): one scan of the two nodes' states, as now, emitting
  runs instead of single ambiguity entries.
- `updateBranchSites` (states changed at given sites only): substitutions as now; the ambiguity
  intervals are only rebuilt if the presence (ambiguous or not) changed at any of the sites, by splicing
  those sites into the interval list: O(runs + changed sites), not O(events).
- `MATInitializer` (parsimony) emits intervals.

### 4. Node summaries instead of event loops in the likelihood

`getComposition` and `getExitRate` recompute a node's composition (summed weights of its sites in each
state) and exit rate λ_n from its parent's and the events on its branch, whenever the node's states
change. Instead, the store keeps for each node, alongside its states:

- `composition[a]`: the summed pattern weights of the node's sites in state a;
- `rateComposition[a]`: the summed weight × site rate of those sites, so that
  λ_n = Σ_a rateComposition[a] · q_a, independent of the rate matrix.

They are updated in O(1) at every state write (through a new `setState(n, p, s)`, replacing direct
writes to the arrays returned by `getWritableStates`), journaled with the states (so restore puts them
back), and computed in full when a node's states are set from scratch (initialisation, rebuild). A
change to a per-site rate updates `rateComposition` at every node, O(nodes) per rate change.

The likelihood then gets λ_n and the composition in O(1) per node, with no loop over events. (This part
does not depend on the intervals and could be done first.)

### 5. The resampler

- Candidate sites: substitutions from the substitution lists, as now. Ambiguity-only sites come as
  *segments*: sweeping the interval end points of the branches around the block gives the maximal ranges
  of patterns over which every node around the block has the same presence (ambiguous or not), in
  O(runs around the block).
- Grouping: within a segment the presence pattern is fixed, so a site's group depends only on its state
  a. The states are counted per segment with one scan of the reference node's state bytes over the
  segment (a byte per site, no per-site conditional, hashing or random draw). Each (presence pattern,
  state) group then contributes count × log P(unchanged), as now.
- Drawing changes: the number of changed sites in a group is drawn by thinning as now; a site to change
  is located by scanning its segment for the k-th site in state a only when a change is drawn (rare).
- Ambiguity flags: after a move only the regrafted parent (and nodes on the path) can change between
  ambiguous and not, and only where its new children's ambiguity intervals intersect differently from
  before, so the writes are found by interval operations, O(runs), instead of rewriting every
  ambiguity-only site's block states.

### 6. Expected effect

The per-site ambiguity work in the profile (about 30% of run time) becomes per-run or per-change work.
I would expect roughly a 20–25% reduction in run time on big-3, more on data with longer runs of N and
less on data with few ambiguities. Memory for the events falls by about two orders of magnitude on
big-3 (355,579 entries to 2,913 intervals).

## Related: cheaper conditionals (not part of this design)

Building the conditionals for each (presence pattern, state) group, before and after each move (about
80 per move), is now the largest single cost (36%): each does a filtering pass in which every node's
partials go through `Math.log`, and `UniformConditional` normalises a row in log space per node. For
groups whose sites keep their state (the vast majority), only P(unchanged) is needed, which can be
computed in linear space with one logarithm per group, and the four states sharing a presence pattern
can share one filtering pass. This is independent of the intervals and may be worth doing first.

## Validation

The existing tests cover the parts that change:

- `testSparseResamplingMatchesDenseConditional`: the Hastings ratio of 300 random moves equals that of
  the dense block conditional, and the MAT likelihood equals the direct path density;
- `testSparseResamplingDrawsUniformSitesExactly`: sites with no events and ambiguity-only sites change
  with their exact probabilities, and ambiguity flags are right after each move;
- the enumeration, consistency and parsimony tests of `MATDataLikelihoodDelegateTest`, and
  `MATBranchHistoryModelTest`.

Plus new unit tests of the interval operations (building, splicing, intersection), and a check that the
store's node summaries always equal those computed from the states. Logs will not be bit-identical to
the current version (the exit rates are summed in a different order, and the resampler draws in a
different order), so the end-to-end check is the chain comparison with the states summed out
(`MATValidation`) and a long big-3 run against the Felsenstein baseline.

## Questions

1. Do the node summaries (section 4) first, as a separate step? They remove the likelihood's event loops
   without intervals and are simpler.
2. Do the cheaper conditionals ("Related") before the intervals? They target a larger share (36%) than
   the intervals (about 30%).
