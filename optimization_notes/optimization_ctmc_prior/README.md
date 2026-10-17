# optimization_ctmc_prior

Branch: `optimization_ctmc_prior` (from `main`)

Only `src/dr/evomodel/tree/CTMCScalePrior.java` changes, plus a new test,
`src/test/dr/evomodel/tree/CTMCScalePriorTest.java`.

## Why

The CTMC scale reference prior (`<ctmcScalePrior>`, used for the clock rate in most BEAUti analyses)
needs the total length of the tree. It cached the tree length, but discarded it whenever anything
changed: the tree, the scale parameter (the clock rate) or the substitution model. It then
recalculated it with `TreeUtils.getTreeLength`, which goes over every node, calls `getBranchLength`
and reads two node heights per branch. With 1,610 tips that costs about 11 µs per evaluation. In a
profile of the coalescent-only benchmark below, this was 9.6% of the run time on `main`.

## Change

- The tree length is no longer discarded when the scale parameter or the substitution model changes,
  since neither affects it (the log density is still recalculated).
- For each tree that is a `DefaultTreeModel`, the prior keeps the length of every branch. It updates
  only the branches that tree changed events say may have changed, following the same rule as
  `TreeDataLikelihood`:
  - a node event (a node's height changed, or it gained or lost a child) marks that node's branch and
    its children's branches;
  - a whole-tree event, or any other event, marks every branch;
  - an event whose object is a node `Parameter` is ignored. `AbstractModel.variableChangedEvent`
    passes on every parameter change as well as the node event that the tree model fires for it,
    so the node event already covers the change.
- The tree length is then the sum of the branch lengths in node order. That is the same sum, in the
  same order, as `TreeUtils.getTreeLength`, so the value is identical to the last bit.
- Branch lengths updated since the last store are journaled, and `restoreState` puts them back. At
  `storeState`, any branches marked but not yet updated are brought up to date first, so that the
  stored lengths match the stored tree.
- Subtree lengths (the `taxa` option) and trees other than `DefaultTreeModel` are still recalculated in
  full, as before.

## Testing

`CTMCScalePriorTest` (JUnit) uses a 40-tip coalescent tree (with some tips sharing a date) and a CTMC
scale prior on a rate. It runs 20,000 store/accept-or-restore cycles of SubtreeLeap, fixed-height SPR,
uniform node height, root height scale, up/down (all internal node heights against the rate) and rate
scale moves. A quarter of the cycles make 2–4 moves before accepting or restoring, the prior is
sometimes evaluated only after several moves, and it is occasionally marked dirty. After every move
and every restore, the prior's log density must equal, exactly, the value calculated from
`TreeUtils.getTreeLength`.

The test passes on `main` and on this branch. Two negative controls fail it: a version that does not
mark the children's branches on a node event, and a version that does not put the journaled lengths
back on restore.

MCMC logs with the same seed are identical (apart from the header lines) to `main`'s, for the test XML
here (10,000,000 states) and for `benchmark-big-2-min` (50,000 states, parameter and tree logs).

## Performance

Machine: Apple M5 Max, OpenJDK 26.0.1, one run at a time, the two runs back to back.

### Coalescent-only benchmark with a CTMC scale prior (`prior-only-ctmc-1610.xml`)

This file is the 1,610 dated EBOV tips from `benchmark-big-2`, with the tree sampled from a
constant-size coalescent and a clock rate under the CTMC scale reference prior, with no sequence data
or tree likelihood. The operators are those of the other optimization notes' `prior-only-1610.xml`
(subtreeLeap, fixedHeightSPR, uniform on node heights, up/down on heights and popSize, popSize scale),
plus a scale on the rate and up/down on all node heights against the rate, as in `benchmark-big-2`.
10,000,000 states, seed 1234.

| CTMCScalePrior | Run time | min / million states | Speed-up |
|---|---|---|---|
| `main` | 20.76 min | 2.08 | – |
| incremental tree length | 19.30 min | 1.93 | **7.0% less time** (1.08×) |

The final log lines are identical (joint −34.5964 at state 10,000,000).

In profiles (2,000,000 states), the CTMC prior's share of the run time fell from 9.6% to 1.7%. Most of
what remains (1.4%) is the sum over the cached branch lengths.

To run:

```
beast -seed 1234 -overwrite prior-only-ctmc-1610.xml
```

### With the tree likelihood (`benchmark-big-2-min`, 50,000 states, BEAGLE)

The run took 44.2 s, against 56.2 s for `main` in an earlier session. The logs are identical. These
runs are short and BEAGLE multi-threaded timings vary by about ±10%, so this does not measure the
change.

## Possible further step

For a binary tree the length can also be written in terms of node heights only:
L = Σ(internal node heights) + root height − Σ(tip heights). A height change of Δ changes the length
by Δ for an internal node, 2Δ for the root and −Δ for a tip, and a topology change with fixed heights
does not change it at all. Updating a running total this way would remove the remaining O(n) sum. But
it would give a different rounding from `TreeUtils.getTreeLength`, so logs would no longer be
bit-identical, and it needs a guard for non-binary trees. It is not worth it for the remaining 1.4%.
