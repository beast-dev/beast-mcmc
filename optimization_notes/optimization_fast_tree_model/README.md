# optimization_fast_tree_model

Branch: `optimization_fast_tree_model`, from `main`, with `optimization_ctmc_prior` merged in (so it
must be merged after that branch).

New: `src/dr/evomodel/tree/FastTreeModel.java` and `src/test/dr/evomodel/tree/FastTreeModelTest.java`.

Changed:

- `src/dr/evomodelxml/tree/TreeModelParser.java`: the `implementation` attribute; returns `TreeModel`.
- `src/dr/evomodel/tree/CTMCScalePrior.java`: the incremental tree length (from
  `optimization_ctmc_prior`) is used with `FastTreeModel` as well as `DefaultTreeModel`;
  `CTMCScalePriorTest` runs with both.
- `src/dr/evomodel/operators/SubtreeSlideOperator.java` and its parser: the operator takes any
  `TreeModel`. It no longer swaps node rates or node traits (on root changes, or the optional random
  swaps); `swapInRandomRate="true"` or `swapInRandomTrait="true"` gives a parse error. With no node
  rates or traits this changes nothing; the few legacy XMLs that combine `subtreeSlide` with node
  traits (e.g. the RacRABV continuous phylogeography examples) will give different (still valid)
  results.
- `src/dr/evomodel/birthdeath/EfficientBirthDeathLikelihood.java`: accepts any `TreeModel` (it only
  checked for a `DefaultTreeModel`).

`DefaultTreeModel` is unchanged.

The design, the decisions taken, and the components that still require a `DefaultTreeModel` are in
[DESIGN.md](DESIGN.md).

## Using it

On `operator_experiment` the fast implementation is the default: a plain `<treeModel>` creates a
`FastTreeModel`. `implementation="old"` selects `DefaultTreeModel`:

```xml
<treeModel id="treeModel">                      <!-- FastTreeModel -->
<treeModel id="treeModel" implementation="old"> <!-- DefaultTreeModel -->
```

(On `main`, where this branch's work was first merged, the default is still `DefaultTreeModel` and the
fast one is selected with `implementation="fast"`.)

The child elements are the same for both (`rootHeight`, `nodeHeights`, `leafHeight`, `leafHeights`), and
references (`<treeModel idref="treeModel"/>`) are unchanged. Node rates and node traits are only
available in the old implementation: a `<treeModel>` with `nodeRates`, `nodeTraits` or `leafTrait` and no
`implementation` attribute selects it automatically, and with `implementation="fast"` they give a parse
error. The few older components that use features of `DefaultTreeModel` (the rateExchange, treeBitMove,
treeBitRandomWalk and RLC narrow exchange operators, the MRCA trait statistic, the inhibition assay
likelihood, the independent coalescent samplers and the star tree model's shared root) give an error
asking for `implementation="old"` if given a fast tree model. The log reports which implementation is
used ("Creating the tree model, 'treeModel' (fast implementation)").

## What it does

`FastTreeModel` behaves exactly as `DefaultTreeModel` (the same node numbering, height parameters and
their ids, bounds and events), so an analysis gives the same output with either. But:

- the tree is held in flat arrays (heights, parent and children, root) and stored copy-on-write, so
  storing, restoring and accepting cost in proportion to the number of nodes changed;
- the node height parameters are views on the arrays, one per `DefaultTreeModel` height parameter
  object ("slot", swapped on root changes as `DefaultTreeModel` swaps its parameter objects); they hold
  no values and are not variables of the model, so storing the model graph does not visit them;
- the compound height parameters (`nodeHeights`) find a changed slot's index directly.

## Testing

`FastTreeModelTest` (JUnit):

- `testLockStepWithDefaultTreeModel` builds a `DefaultTreeModel` and a `FastTreeModel` from the same
  11-tip tree (with tips at different heights), with root, internal, all-internal and two leaf height
  parameters (one with an extra bound, as the parser adds for an uncertain tip date). It then makes
  20,000 store / move(s) / accept-or-restore cycles, applying the same move with the same random numbers
  to both: subtree leap, fixed-height SPR, narrow and wide exchange, Wilson–Balding, uniform on internal
  heights, scale on the root height, up/down on all internal heights, a random walk on one tip height
  and uniform on the other. After every move and every restore it checks that the trees (structure,
  child order, heights, root, taxa) and the values, bounds and names of every height parameter are
  identical, and at the end that the full sequence of events (tree events and parameter events, as
  listeners see them) is identical. There are over 100 root changes and over 100,000 events. It passes
  with and without assertions enabled (with them, the tree validity check at the end of each edit runs
  too). Two negative controls fail it: a version that does not swap slots on a root change, and one
  that does not restore the slots.
- `testNodeRatesAndTraitsUnsupported` checks that node rate and trait methods throw.

MCMC logs with the same seed are identical (apart from the header lines) between
`implementation="fast"` and the default, with assertions enabled, for:

| XML | States | Default | Fast |
|---|---|---|---|
| coalescent-only with CTMC prior (`prior-only-ctmc-1610.xml`) | 200,000 | 25.3 s | 11.6 s |
| the same plus narrow and wide exchange, Wilson–Balding, root height scale, a sampled tip date and a tip with an uncertain date | 500,000 | 1.36 min | 34.7 s |
| the same plus a birth-death tree prior (`Gernhard08BirthDeathModel`) | 500,000 | 1.52 min | 42.9 s |
| `benchmark-big-2-min` (BEAGLE tree likelihood; parameter and tree logs) | 50,000 | – | 39.9 s |
| the operators XML plus `subtreeSlide` (default run with `main`'s operator) | 500,000 | 1.32 min | 35.7 s |
| the birth-death XML plus `birthDeathLikelihood` (an episodic birth-death-sampling model; default run with `main`'s likelihood) | 500,000 | 2.05 min | 1.18 min |

and for the 10,000,000-state runs below.

## Performance

Machine: Apple M5 Max, OpenJDK 26.0.1, one run at a time.

### Coalescent-only benchmark with a CTMC scale prior (`prior-only-ctmc-1610.xml`)

The 1,610 dated EBOV tips from `benchmark-big-2`, with the tree sampled from a constant-size coalescent
and a clock rate under the CTMC scale reference prior (no sequence data or tree likelihood).
10,000,000 states, seed 1234. "Combined" is `main` with the other optimization branches merged
(`optimization_fast_intervals`, `optimization_tree_store` (with `optimization_tree_restore`),
`optimization_ctmc_prior`, `optimization_compound_parameter`, `optimization_parameter_store`).

| Build | Tree model | Run time | min / million states | vs `main` default |
|---|---|---|---|---|
| `main` | default | 20.96 min | 2.10 | – |
| `main` + FastTreeModel (without the CTMC change) | fast | 9.73 min | 0.97 | 54% less time (2.2×) |
| this branch | fast | 8.82 min | 0.88 | **58% less time** (2.4×) |
| combined | default | 6.95 min | 0.70 | 67% less time (3.0×) |
| combined + this branch | fast | 4.14 min | 0.41 | **80% less time** (5.1×) |

The two rows for FastTreeModel on `main` were run back to back (run times on this machine vary by up
to about 5% between sessions). All final log lines are identical (joint −34.5964 at state 10,000,000).

On its own, FastTreeModel gives about the same saving on `main` as the three tree and parameter store
branches together, without changing `DefaultTreeModel` or `Parameter`.

### Where the time goes now (profiles of 2,000,000 states)

With the combined branches and FastTreeModel (before the CTMC prior used its incremental tree length
with FastTreeModel), storing the model graph was 13% of the run time (restore and accept under 1%), the
CTMC prior's full tree length 18%, the coalescent likelihood 52% and the operators 7%. Enabling the
CTMC prior's incremental tree length removed most of that 18% (5.24 → 4.14 minutes).

On this branch alone (with `main`'s FastIntervals), sorting the coalescent intervals takes over half of
the run time (`optimization_fast_intervals` addresses this) and the CTMC prior 3%.

To run:

```
beast -seed 1234 -overwrite prior-only-ctmc-1610-fast.xml
```

(`prior-only-ctmc-1610.xml` is the same XML with the default tree model.)
