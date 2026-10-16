# optimization_fast_tree_model

Branch: `optimization_fast_tree_model` (from `main`)

New: `src/dr/evomodel/tree/FastTreeModel.java` and `src/test/dr/evomodel/tree/FastTreeModelTest.java`.
Changed: `src/dr/evomodelxml/tree/TreeModelParser.java` (the `implementation` attribute; returns
`TreeModel`). `DefaultTreeModel` is unchanged.

The design, the decisions taken, and the components that still require a `DefaultTreeModel` are in
[DESIGN.md](DESIGN.md).

## Using it

```xml
<treeModel id="treeModel" implementation="fast">
    ...
</treeModel>
```

The child elements are the same as for the default implementation (`rootHeight`, `nodeHeights`,
`leafHeight`, `leafHeights`), and references (`<treeModel idref="treeModel"/>`) are unchanged. Node
rates and node traits are not supported: `nodeRates`, `nodeTraits` and `leafTrait` give a parse error,
and analyses that need them must use the default implementation.

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

and for the 10,000,000-state runs below.

## Performance

Machine: Apple M5 Max, OpenJDK 26.0.1, one run at a time, all four runs back to back.

### Coalescent-only benchmark with a CTMC scale prior (`prior-only-ctmc-1610.xml`)

The 1,610 dated EBOV tips from `benchmark-big-2`, with the tree sampled from a constant-size coalescent
and a clock rate under the CTMC scale reference prior (no sequence data or tree likelihood).
10,000,000 states, seed 1234. "Combined" is `main` with the other optimization branches merged
(`optimization_fast_intervals`, `optimization_tree_store` (with `optimization_tree_restore`),
`optimization_ctmc_prior`, `optimization_compound_parameter`, `optimization_parameter_store`).

| Build | Tree model | Run time | min / million states | vs `main` default |
|---|---|---|---|---|
| `main` | default | 20.96 min | 2.10 | – |
| `main` | fast | 9.01 min | 0.90 | **57% less time** (2.3×) |
| combined | default | 7.08 min | 0.71 | 66% less time (3.0×) |
| combined | fast | 5.24 min | 0.52 | **75% less time** (4.0×) |

All four final log lines are identical (joint −34.5964 at state 10,000,000).

On its own, FastTreeModel gives about the same saving on `main` as the three tree and parameter store
branches together, without changing `DefaultTreeModel` or `Parameter`.

### Where the time goes now (combined, fast; profile of 2,000,000 states)

| | Share of run time |
|---|---|
| storing the model graph (`storeModelState`) | 13% (restore and accept < 1%) |
| CTMC scale prior: full tree length (`TreeUtils.getTreeLength`) | 18% |
| coalescent likelihood (TreeIntervals adding every node's height, FastIntervals, the likelihood itself) | 52% |
| operators | 7% |

The CTMC prior recalculates the whole tree length here because its incremental cache (branch
`optimization_ctmc_prior`) is only enabled for `DefaultTreeModel`. FastTreeModel fires the same events,
so enabling it for FastTreeModel too would remove most of that 18%.

To run:

```
beast -seed 1234 -overwrite prior-only-ctmc-1610.xml
```

(`prior-only-ctmc-1610-fast.xml` is the same XML with `implementation="fast"`.)
