# The MAT resamples its history when evaluated, and reports the Hastings term itself

Branch `operator_experiment`, October 2026.

## What changed

Before, the mutation-history resampling that accompanies a subtree move was a `SubtreeMoveHook` attached in
XML to each subtree operator; the operator called it before and after its tree edit and added its log
Hastings ratio to its own. Now the resampler is attached to the likelihood, and:

- `MATDataLikelihoodDelegate` finds the subtree move when it is next evaluated after a topology change,
  from the nodes whose parent in the tree differs from their parent in the store (found among the nodes
  with a tree event since the step began, so O(touched), not O(n)); runs the resampler's old phase through
  a `TopologyView` of the tree before the move (the store's parents, the children of the few nodes whose
  children changed, and the heights of the last accepted state, which the delegate keeps) and its new phase
  through the tree; and accumulates the log Hastings ratio.
- `MarkovChain` collects that ratio after evaluating the posterior through the new
  `dr.inference.model.HastingsRatioContributor` interface (implemented by `TreeDataLikelihood`, delegating
  to its delegate) and adds it to the operator's ratio before the acceptance test. The term is therefore
  not tempered under MC3, and nothing about operators changed: `MCMCOperator.operate()` still returns a
  double.
- The `SubtreeMoveHook` interface, the hook plumbing in `AbstractTreeOperator` and the calls in
  `SubtreeLeapOperator`, `SubtreeJumpOperator` and `FixedHeightSubtreePruneRegraftOperator` are gone. A
  `<mutationHistoryResampler>` child on an operator is a parse error that says where to put it.

Any operator that makes a single subtree prune-and-regraft (with or without a height change, at the root
or not) now works with the MAT: subtree leap, FHSPR, subtree jump, Wilson-Balding, subtree slide, narrow
exchange, NNI, FNPR. A change that is not a single SPR (wide exchange, `GibbsSubtreeSwap`, a whole-tree
replacement, or two topology changes in one step) stops the run with a clear message when a resampler is
attached. Outside an MCMC step (a tree loaded from a checkpoint) the MAT is rebuilt by parsimony as before.

### The canonical description of a move

When the moved parent P lands on the branch below its own sibling, or above its own parent, the same change
of topology is also a legal move of a different subtree, resampling a different block (P at all sites plus
the sibling at the hot sites, or the other way round). The delegate must choose the same description for a
move and its reverse, or detailed balance fails silently. It chooses the description with the smallest
moved-parent number, a function of the pair of topologies alone. Under the old hook the operator fixed the
description, so this did not arise. `testLazyMoveDescriptionIsReversible` checks it over random moves,
including those shapes and moves at the root.

### Limits

- One resampler per MAT likelihood (one block type); before, each operator could have its own.
- A strict clock (or none): the branch rates of the old topology are read from the current tree by node,
  which a relaxed clock's node-to-rate mapping does not allow after a root change. `setResampler` rejects
  other branch rate models.
- Operators that evaluate the likelihood for several candidate topologies in one `operate()`
  (`GibbsPruneAndRegraft`, `GibbsSubtreeSwap`) are unsupported with a resampler (the second topology change
  in a step throws).

## XML

```xml
<matDataLikelihood id="treeLikelihood"> ... </matDataLikelihood>
<mutationHistoryResampler id="matResampler" block="local">
    <matDataLikelihood idref="treeLikelihood"/>
</mutationHistoryResampler>
<operators>
    <subtreeLeap size="1.0" weight="165"><treeModel idref="treeModel"/></subtreeLeap>     <!-- no resampler child -->
    <wilsonBalding weight="16"><treeModel idref="treeModel"/></wilsonBalding>              <!-- now works too -->
    <mutationHistoryGibbsOperator weight="16" size="3">
        <mutationHistoryResampler idref="matResampler"/>                                   <!-- unchanged -->
    </mutationHistoryGibbsOperator>
</operators>
```

The 21 `-mat` XMLs of the beast-benchmark repository were updated (the resampler children removed from
`subtreeLeap` and `fixedHeightSubtreePruneRegraft`); the XML-generating scripts under `optimization_notes`
needed no change.

## Checks

- `ant junit` classes `MATDataLikelihoodDelegateTest` (25 tests: the exact-conditional and Hastings-ratio
  identities now go through the lazy path; new: reversible move description, unsupported change refused in a
  step and rebuilt outside one, height-only moves contribute nothing), `MATBranchHistoryModelTest` (6) and
  the new `MarkovChainTest` (a stub contributor's ratio is added to the operator's, untempered, taken
  once). All pass with `-ea`.
- small-2 MAT, 300k states with `fullEvaluation="5000"` and `-ea`: no evaluation or restore errors, no
  assertion (the delegate asserts at every store that its kept heights match the tree), acceptance rates as
  before. The same with `wilsonBalding` and `narrowExchange` added (accepted 3% and 37%); with
  `wideExchange` the run stops at the first step with the message above.
- Speed, big-3 local block, 1M states, one core: `subtreeLeap` step 62.3 µs (62.0 before); the resampling
  moved from the propose column (29.8 → 2.7 µs) to the evaluate column (28.7 → 55.8 µs). 1.18 min per
  million states, as before.
- Posterior equivalence: see below.

## Posterior validation (small-2, 50M states, seeds 11 22 33 44)

Method `mat-lazy-gibbs` of `optimization_notes/mat_local_spr/benchmark/make_validation_xmls.py` (local
block, `subtreeLeap` + `fixedHeightSubtreePruneRegraft` without resampler children, the Gibbs operator as
before) against the BEAGLE `baseline`, four seeds each, 10% burn-in, analysed with `analyse.py`. Full
tables in [validation-results.md](validation-results.md).

- Every within-run, between-replicate and method-versus-reference z-score is below 3 in magnitude. Pooled
  over seeds, the lazy MAT differs from BEAGLE by |z| <= 0.6 on rootHeight, treeLength, popSize,
  clock.rate and the coalescent.
- Posterior means agree with the earlier hook-based `mat-local-gibbs` runs
  (`optimization_notes/mat_local_spr/benchmark/validation-results.md`) and with BEAGLE:

| Parameter | BEAGLE (4 seeds) | lazy MAT (4 seeds) | hook-based MAT (earlier) |
|---|---|---|---|
| rootHeight | 6.202-6.214 | 6.189-6.215 | 6.176-6.220 |
| treeLength | 62.48-62.68 | 62.33-62.70 | 62.16-62.69 |
| constant.popSize | 3.696-3.730 | 3.689-3.728 | 3.687-3.728 |
| clock.rate (x1e-3) | 3.718-3.739 | 3.717-3.741 | 3.716-3.754 |
| treeLikelihood | - | -5785 | -5785 |
| mat.mutations | - | 391.88-392.02 | 391.91-391.97 |
| excess mutations | - | 0.313-0.345 | 0.310-0.343 |

- ESS per million states for the lazy MAT (rootHeight 1117-1483, treeLength 535-842) matches the hook-based runs
  (1079-1461 and 519-720); the proposal is unchanged, so this is expected. Wall time 0.29-0.30 min per million
  states, the same as the earlier hook-based figure (0.29-0.32 for the two unloaded seeds). The BEAGLE
  baseline took 1.05-1.10 min per million in this study (0.43-0.85 earlier); wall times are not comparable
  across studies as the machine load differed.
