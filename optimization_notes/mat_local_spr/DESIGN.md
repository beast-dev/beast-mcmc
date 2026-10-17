# MAT: a local subtree-move hook (design)

Branch: `operator_experiment`. Status: design, for agreement before implementation.

The target is speed (effective samples per hour) on sparse mutation trees of 1,000 to 10,000 tips.

## Why

`MutationHistoryResampler` is attached to subtree moves (`subtreeLeap`, `fixedHeightSubtreePruneRegraft`)
through `SubtreeMoveHook`. Each move does an exact Gibbs block update of the states of a *path block*:

- the regrafted parent P;
- every internal node on the path between the old sibling S and the destination D in the pruned tree.

This covers every site with an event around that block. On big-3 (899 BDBV genomes, integrated hits),
after the optimisations in commits 92f324168 to 1ba9bbaae, the run takes 38 s per 300,000 states
(about 127 µs per state). The resampler is about 65% of that. Per move it handles:

- about 1,150 sites whose only events around the block are ambiguities, in about 27 groups;
- the conditional of each group, built twice (under the old and the new topology);
- a full scan of the patterns of each reparented branch.

All of this grows with the length of the path, so it grows with the leap distance.

Varilly et al. (Nature 2026) describe subtree prune-and-regraft moves whose cost scales with the
mutations near the moved subtree instead:

- the subtree carries its own mutations;
- only the sites at which the subtree is the only data along part of its path to the rest of the tree
  are re-sampled along that path;
- missing data are held as intervals of sites, so finding those sites costs per run of missing sites,
  not per site.

This design brings that locality into the MAT package. The subtree moves and the `SubtreeMoveHook`
interface stay as they are.

## Interface (unchanged)

```java
void beforeSubtreeMove(NodeRef node, NodeRef sibling, NodeRef destination); // old topology
double afterSubtreeMove();                                                  // new topology; log Hastings ratio
```

The operator chooses the move and adds the hook's log Hastings ratio to its own. Nothing outside the MAT
package changes.

Notation (X = `node`):

| | Before the move | After the move |
|---|---|---|
| P, the parent of X | between S and G = parent(P) | between D and D' = parent(D) before the move |
| S | child of P | child of G |
| D | child of D' | child of P |

## The local block

### Presence is deterministic

The MAT keeps the invariant that a node has a state at a site (is *present*) if and only if some tip below
it has data there. Moving X changes presence only at ancestors whose only data at a site comes from X:

- on the old side, from G upwards;
- on the new side, from D' upwards.

These ancestors lie on the path between the old and new positions, below their most recent common
ancestor; above that, X is below both before and after. The current path block contains them all, which
is why it is that large.

### Hot sites

A site l is *hot at level k* on the new side when:

- X has data at l;
- D', and the k − 1 ancestors above it, have no data at l before the move.

Those k nodes become present only through X. At a hot site, the states of the chain of newly present
ancestors, down to P, have to be drawn: between the first ancestor T that already had data (whose state
is fixed) and X (whose state is fixed).

The old side is the mirror image. At sites where G and the ancestors above it have no data other than X,
those nodes become absent, which is deterministic.

The hot sites at each level are the intersection of the missing runs of the siblings met going up the
path (S, then the sibling at each level). The resampler already caches each node's ambiguous runs by
state version (d99267eda). So this is an intersection of interval lists that stops as soon as it is empty,
usually within a level or two: O(runs on the siblings), not O(sites).

### What is resampled

The block is:

- P's states at all variable sites;
- at each hot site, the states of the chain of newly present ancestors (the new side) or of the ancestors
  that lose their data (the old side).

No other state changes. In particular, the internal nodes on the path keep their states wherever they
have data from elsewhere.

- **Non-hot sites:** P's state is drawn from its full conditional given its three neighbours (D', D and X
  after the move; G, S and X before). This is the existing `Conditional` with a block of one node.
- **Hot sites:** the chain's states are drawn from their full conditional given T and X: forward filtering
  and backward sampling on a chain, with the existing per-branch matrices.

### Hastings ratio

The reverse move (X back from D to S) resamples the same kinds of variables:

- P's states at all sites;
- at the old side's hot sites, the chain that becomes present again, which is the chain the forward move
  made absent.

Each direction draws its block from the exact full conditional given everything else, under its own
topology. So the ratio is that of the current hook:

```
log q(old block states | old topology) − log q(new block states | new topology)
```

`beforeSubtreeMove` scores the old block, and `afterSubtreeMove` draws the new one and returns the
difference. The presence changes are functions of the topology and add nothing to the ratio.

## Sites and cost

The candidate sites are:

- the events on the branches incident to P before and after the move: the branches above X, S, P (old),
  D and P (new);
- the hot sites on both sides.

The other variable sites are uniform around P, and are handled by counts and thinning as now (`sampleUniformSites`):

- their state a is read from X, whose states do not change;
- with three neighbours, there are at most four of these conditionals, one per state.

Sites whose only events on those branches are ambiguities are grouped by:

- which of P's three neighbours is present (at most 8 patterns);
- the state.

The groups are found by the run sweep added in d99267eda, over three nodes instead of a whole path.

The branch events can also be updated without full scans:

- **S's new branch (S → G)** combines the events of S → P and P → G site by site.
- **D's branch** splits into D' → P and P → D. They differ only at the candidate sites where P's new
  state differs from D's or D''s.
- **X's branch** changes only at candidate sites.

Each update is O(events), which removes the full `updateBranch` pattern scans (about 6% of the run time
now).

So a move costs O(events on five branches + runs on the siblings met + changes). That doesn't depend on
the leap distance. What remains per site is the ambiguity-only sites on those few branches:

- a tip can have hundreds of missing sites in a handful of runs;
- if that still shows in the profile, those sites are counted per run (segment counting, or the
  intervals of `mat_ambiguity_intervals/DESIGN.md`) instead of per site.

## Branch history models

- **Integrated and single hits:** as above. These are the first target.
- **Explicit hits:** the ancestral states are drawn as above, then the hits on the branches incident to P
  and on the hot chains are drawn given their end states, as the current hook does for its block. The
  probabilities of the old and new hits enter the Hastings ratio.
- **Timed hits:** inserting P on D's branch splits that branch's hits by time, as in Varilly et al.
  (2026). D's hits do not need to be resampled.

## Mixing

Path nodes keep their states, so a long move puts all the differences between X and its new
neighbourhood onto the branch from P to X:

- acceptance falls with distance;
- `subtreeLeap` adapts its kernel towards its target acceptance (`targetAcceptance`, default 0.234), so
  the leaps will be shorter but much cheaper.

Whether 0.234 is still a good target per hour (rather than per iteration) is an experiment, below.

The chain stays irreducible:

- every internal node is P for some move;
- under the branch history models each single-node conditional gives positive probability to every
  state.

But a mutation can no longer move across several nodes in one step. So the path hook remains available
for wide moves.

## Configuration

The block is a property of each hook element, so different operators can use different hooks on the
same MAT likelihood:

```xml
<mutationHistoryResampler id="localResampler" siteSelection="all" block="local">
    <matDataLikelihood idref="treeLikelihood"/>
</mutationHistoryResampler>
<mutationHistoryResampler id="pathResampler" siteSelection="all" block="path">
    <matDataLikelihood idref="treeLikelihood"/>
</mutationHistoryResampler>

<subtreeLeap size="1.0" weight="899">
    <treeModel idref="treeModel"/>
    <mutationHistoryResampler idref="localResampler"/>
</subtreeLeap>
<subtreeLeap size="5.0" weight="90" targetAcceptance="0.1">
    <treeModel idref="treeModel"/>
    <mutationHistoryResampler idref="pathResampler"/>
</subtreeLeap>
<fixedHeightSubtreePruneRegraft weight="89.9">
    <treeModel idref="treeModel"/>
    <mutationHistoryResampler idref="localResampler"/>
</fixedHeightSubtreePruneRegraft>
```

- `block="path"` is the current behaviour and stays the default until the experiments below are done.
- Two `subtreeLeap`s with the same hook and target acceptance would adapt to the same kernel. A wide
  component needs a fixed kernel or a different `targetAcceptance`.
- `fixedHeightSubtreePruneRegraft` is mainly a burn-in aid (wide moves with the heights, and so the
  coalescent density, fixed). Which hook suits it is an experiment too.

## Validation

- **Hastings ratio:** on small trees with ambiguous sites, compare the log Hastings ratio of random moves
  with that of the dense full conditional of the local block, found by enumeration (generalising
  `testSparseResamplingMatchesDenseConditional`).
- **Exact draws:** P's states and the hot-chain states change with their exact probabilities, and
  presence is right after every move (generalising `testSparseResamplingDrawsUniformSitesExactly`).
- **Invariant:** after every move, a node is present at a site if and only if a tip below it has data
  there, and each branch's events equal those rebuilt from the states.
- **Chains:** with the states sampled, chains agree with the states summed out (`MATValidation`). The MCMC
  full-evaluation check (`-Dmcmc.evaluation.count`, threshold 1e-6) passes.
- **Bias:** a long big-3 run agrees with the Felsenstein baseline.

## Experiments

1. **Sampling:** ESS per hour (not per iteration) of the root height, the tree length, the coalescent
   parameters and split frequencies.
   - Hooks: `local` and `path`, each alone and in the mixture above.
   - Target acceptances: 0.1, 0.234 and 0.4.
   - Data: big-3 and a sparse data set of about 10,000 tips.
2. **Burn-in:** time to stationarity from the same starting tree with `fixedHeightSubtreePruneRegraft`
   using each hook.

## Implementation steps

1. Hot sites by intersecting cached ambiguous runs, with unit tests against a direct computation.
2. The local block sampler: P's conditional (the existing `Conditional` with one node), the hot chains, and
   the old-side scoring.
3. Branch-event updates in O(events) for the reparented branches.
4. The `block` attribute, and tests 1–4 of the validation.
5. Experiments.
6. Per-run counting of ambiguity-only sites, if the profile still calls for it.

## Questions

1. Keep `block="path"` as the default until the experiments, or switch the default once validated?
2. Hot sites are drawn per site from a conditional shared by a group: (level, state of T, state of X). Most
   of them keep their state, so they would be drawn by thinning as the uniform sites are. Is it worth
   doing that from the start, or only if they show in the profile?
3. Should a hook ever influence the choice of destination? Varilly et al. (2026) choose destinations near
   the moved subtree. Here the operator's distance kernel does that, which keeps the operators unaware of
   the MAT.
