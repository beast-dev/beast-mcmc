# FastTreeModel: design

Branch: `optimization_fast_tree_model` (from `main`). Status: agreed and implemented (see Decisions at the end).

## Aim

A new tree model, `dr.evomodel.tree.FastTreeModel`, for standard analyses on large trees. It is a
drop-in alternative to `DefaultTreeModel` (selected in XML with `<fastTreeModel>` in place of
`<treeModel>`) that:

- gives **bit-identical** MCMC logs to `DefaultTreeModel` for the same XML and seed;
- makes storing, restoring and accepting the tree cost O(nodes changed), not O(nodes), with no
  per-node objects visited when the model graph is stored;
- leaves `DefaultTreeModel` unchanged, so nothing that depends on it can break.

Node rates and node traits are not supported: every `TreeModel` method for them throws
`UnsupportedOperationException`, and the parser rejects the elements that create them. Analyses that
need them (the autocorrelated relaxed clock, sampled ancestral traits, tip traits for continuous-trait
models) keep using `DefaultTreeModel`.

## Why DefaultTreeModel is slow on large trees

Profiles of the coalescent-only benchmarks (1,610 tips, 3,219 nodes), on `main` and with the other
optimization branches combined, show that the remaining tree cost comes from having one real
`Parameter` object per node:

- every node height `Parameter` is added to the tree model with `addVariable`, so every store, restore
  and accept of the model graph visits all 3,219 of them (with the other branches combined, about 14%
  of the run time just visiting them, plus part of the 12% in the `AbstractModel` walks);
- each `Parameter` stores its own values (copied eagerly on `main`);
- the tree keeps a map from parameter to node, rebuilt on every restore on `main`;
- the tree keeps a second array of nodes and copies every node's links at every store on `main`;
- on a root change, `swapParameterObjects` exchanges the parameter objects of the old and new root, so
  that the root height parameter is always the root's.

The other branches (`optimization_tree_restore`, `optimization_tree_store`,
`optimization_parameter_store`, `optimization_compound_parameter`) reduce these costs but cannot remove
the per-node visits without changing how every model is stored.

## State

All tree state is held in flat arrays indexed by node number (tips `0 … n−1`, internal nodes
`n … 2n−2`, as in `DefaultTreeModel`):

| Array | Contents |
|---|---|
| `double[] heights` | node heights |
| `int[] parent`, `int[] child0`, `int[] child1` | links (−1 for none) |
| `int root` | the root's node number |
| `int[] nodeSlot`, `int[] slotNode` | which height-parameter slot each node holds, and the inverse (see below) |

Node objects (`NodeRef`) are immutable: a number and, for tips, a taxon. They are created once, so a
`NodeRef` stays valid for the life of the model (unlike `DefaultTreeModel`, whose node objects were
replaced by their stored copies on every restore).

## Store, restore and accept: copy on write

As in `optimization_tree_store`, but for the arrays:

- `storeState` records the root number and clears the journal of changed nodes, O(changed).
- The first change to a node after a store (its height, any of its links, or its slot) saves that
  node's height, links and slot, and adds it to the journal.
- `restoreState` puts back the saved values of the journaled nodes and the root. O(changed).
- `acceptState` does nothing (the journal is cleared at the next store).

The height parameters (below) hold no values of their own and are not added with `addVariable`, so
the model graph walk visits only the tree model itself.

## Height parameters: slots and views

### Slots, for bit-identical behaviour

In `DefaultTreeModel` each node owns a height `Parameter` object, and the compound parameters
(`nodeHeights`) list those objects in node-number order at the time they are created. On a root change
the old and new root exchange parameter objects (keeping their heights), so element *k* of a compound
parameter can refer to a different node afterwards. Operators choose elements by index, so to draw
the same random numbers and move the same nodes as `DefaultTreeModel`, FastTreeModel must reproduce
this mapping exactly.

A **slot** stands for one of `DefaultTreeModel`'s height parameter objects. Slot *s* starts at node
*s*. `slotNode[s]` is the node that slot is now at, and `nodeSlot[i]` the slot at node *i*. On a root
change (`endTreeEdit` with a different root), the slots of the old and new root are exchanged, as
`swapParameterObjects` exchanges the parameter objects. Slot changes are journaled like everything
else.

### Per-slot height parameters

`FastTreeModel` creates one small parameter per slot, a `SlotHeightParameter` (a one-dimensional
`Parameter.Proxy`):

- `getParameterValue(0)` returns `heights[slotNode[s]]`;
- `setParameterValue(0, v)` asks the tree to set the height of `slotNode[s]` (journal, write the
  array, fire events: see below);
- store, restore and accept do nothing (the tree's journal covers them);
- its bounds are those of `DefaultTreeModel`'s `NodeHeightBounds` (the parent's height above, the
  higher child's height below, 0 below a tip, no upper limit for the root), intersected with any bounds
  added later (the parser adds bounds for tips with uncertain dates);
- listeners (priors on `treeModel.rootHeight`, compound parameters, operators) are attached to it in
  the usual way.

These are what `getRootHeightParameter()` (the slot at the root, which by construction stays at the
root) and `getLeafHeightParameter(node)` return. They are ordinary objects, one per node, but unlike
`DefaultTreeModel`'s they are not model variables, so they are never visited at store, restore or
accept, and they hold no values to copy.

### Multi-node height parameters

`createNodeHeightsParameter(rootNode, internalNodes, leafNodes)` returns a `NodeHeightsParameter` (a
multi-dimensional `Parameter.Proxy`) over a list of slots, in the same order as `DefaultTreeModel`
builds its compound parameter (internal nodes in number order, filtered by `rootNode`/`internalNodes`
at creation, then the tips). Element *k* is slot `slots[k]`; it keeps an `int[]` from slot to element
so that a change to a slot is passed on with its index in O(1) (no search, unlike `CompoundParameter`
on `main`). Its bounds are the per-slot bounds. Setting element *k* sets that slot's height.

`<leafHeights>` (a compound of selected tips' leaf heights) is built the same way, over the selected
tips' slots.

## Events

The events must match `DefaultTreeModel`'s, because tree likelihoods, interval models and priors rely
on them (and their order can matter for bit-identical logs, for example when a listener calculates
something as the event arrives).

When a node height changes, `DefaultTreeModel` does this: the node's height `Parameter` sets its value
and fires a parameter-changed event to its listeners, in the order they were added; the tree model is
always the first listener, so it pushes a `TreeChangedEvent(node, parameter, index, height=true)`
(fired immediately, or held until `endTreeEdit` inside an edit), and `AbstractModel` then also fires a
model-changed event with the parameter as the object; then the other listeners of the parameter (for
example compound parameters, which pass the change on with their own index) are notified. (The order
of listeners is to be confirmed against `DefaultTreeModel` in the lock-step test below.)

FastTreeModel reproduces this sequence for the slot of the node that changed: the tree's own
notifications first, then the slot parameter's other listeners, then each multi-node view that
contains the slot, in the order those were created.

Topology changes fire the same events as `DefaultTreeModel`: `addChild` pushes a node event for the
parent, `setRoot` one for the new root, and `beginTreeEdit`/`endTreeEdit` hold and then fire them.
Whole-tree changes (`adoptTreeStructure`) fire a whole-tree event.

## TreeModel API

FastTreeModel extends `TreeModel` and implements the same `Tree`/`MutableTree` methods as
`DefaultTreeModel` from the arrays: node counts, `getNode`, `getRoot`, `getParent`, `getChild`,
`getChildCount`, `isRoot`, `isExternal`, `getNodeHeight`, `getNodeHeightUpper`/`Lower`,
`setNodeHeight`, `setNodeHeightQuietly`, `setRoot`, `addChild`, `removeChild`, `beginTreeEdit`,
`endTreeEdit`, `isTreeValid`, `getNodeTaxon`, taxa, `adoptTreeStructure` (checkpointing),
`isTipDateSampled`.

These throw `UnsupportedOperationException`:

- node rates: `setNodeRate`, `createNodeRatesParameter`; `hasRates` returns false. (`getNodeRate`
  returns 1.0, and `getNodeAttribute`/`getNodeAttributeNames` behave as in a `DefaultTreeModel` without
  node rates, because generic tree copying, such as `new FlexibleTree(tree)` and `SimpleNode`, asks
  every tree for its node rates);
- node traits: `getNodeTrait`, `setNodeTrait`, `getMultivariateNodeTrait`, `setMultivariateTrait`,
  `getNodeTraitParameter`, `getTraitMap`, `swapAllTraits`, `createNodeTraitsParameter`,
  `createNodeTraitsParameterAsMatrix`, `hasNodeTraits` returns false;
- `setNodeTaxon` and `replaceChild` (as in `DefaultTreeModel`, `replaceChild` is unimplemented).

Only binary trees are supported (the parser resolves the starting tree, as `DefaultTreeModel` does).

## XML

```xml
<fastTreeModel id="treeModel">
    <coalescentTree idref="startingTree"/>
    <rootHeight>
        <parameter id="treeModel.rootHeight"/>
    </rootHeight>
    <nodeHeights internalNodes="true">
        <parameter id="treeModel.internalNodeHeights"/>
    </nodeHeights>
    <nodeHeights internalNodes="true" rootNode="true">
        <parameter id="treeModel.allInternalNodeHeights"/>
    </nodeHeights>
    <leafHeight taxon="...">
        <parameter id="age(...)"/>
    </leafHeight>
</fastTreeModel>
```

`FastTreeModelParser` accepts the same attributes (`fixHeights`, `fixTree`) and the same child elements
as `TreeModelParser` for heights: `rootHeight`, `nodeHeights`, `leafHeight`, `leafHeights` (including
the bounds for uncertain tip dates). `nodeRates`, `nodeTraits` and `leafTrait` are rejected with an
error saying to use `<treeModel>`.

References need care. A BEAUti XML refers to the tree model as `<treeModel idref="treeModel"/>` in many
places. `XMLParser` accepts a reference whose element name differs from the referenced element only if
the referring element's parser returns a compatible class, and `TreeModelParser.getReturnType()` returns
`DefaultTreeModel.class`. So with a separate `<fastTreeModel>` element, each `<treeModel idref=…/>`
would give a warning (and an error with `-strict`) unless it were also renamed. Two ways round this:

- change `TreeModelParser.getReturnType()` to `TreeModel.class` (a one-line change, which also matches
  what the parsers that take a tree actually ask for); or
- select the implementation with an attribute on `<treeModel>` instead of a new element (question 1).

Code that casts to `DefaultTreeModel` (about 15 places outside the tree package: some operator parsers,
older continuous-trait and speciation models) will not accept a FastTreeModel. These are the legacy
cases that use node rates or traits, or should be checked individually.

## Validation

1. Unit tests, as for `optimization_tree_store`: random moves (SubtreeLeap, fixed-height SPR, narrow and
   wide exchange, Wilson–Balding, uniform, scale, up/down, tip-date moves), store/accept/restore, several
   moves before a restore and root changes; checking the structure, heights, slot mapping, bounds and
   that every restore reproduces the stored tree exactly.
2. A lock-step test: the same moves applied to a `DefaultTreeModel` and a `FastTreeModel` with the same
   seed, comparing the trees, the parameter values seen through every height parameter and the
   sequence of events after every step.
3. Bit-identical MCMC logs against `<treeModel>` for: the coalescent-only benchmarks, `benchmark-big-2`
   (with the tree likelihood), a birth-death XML, and an XML with tip-date sampling (`leafHeight`
   operators and uncertain dates).
4. Timing against `DefaultTreeModel`, on `main` and with the other optimization branches.

## Out of scope

- Node rates and node traits (use `DefaultTreeModel`).
- Non-binary trees (`StarTreeModel` and similar keep their own classes).
- Changing BEAUti to write `<fastTreeModel>`: a separate decision once FastTreeModel is validated.

## Questions

1. Name of the XML element: `<fastTreeModel>`, or an attribute on `<treeModel>` (e.g.
   `<treeModel implementation="fast">`)? The attribute keeps one element name, so existing
   `<treeModel idref=…/>` references work unchanged (see XML above), and switching an XML is a
   one-attribute edit. `TreeModelParser` would then create either class, and its return type would
   become `TreeModel`.
2. Should the per-slot `SlotHeightParameter`s be created lazily (only when a parser or operator asks
   for a height parameter)? Most analyses ask for the root and all internal heights anyway, so this
   would mainly help analyses that use only tree moves.

## Decisions

1. The implementation is chosen with an attribute: `<treeModel id="treeModel" implementation="fast">`
   (`implementation="default"` or no attribute gives `DefaultTreeModel`). `TreeModelParser` creates
   either class and its return type is now `TreeModel`, so `<treeModel idref=…/>` references work for
   both. The parser may later default to the fast implementation.
2. The per-slot height parameters are created when the tree model is created: most XMLs use the
   up/down operator on all node heights anyway.
3. Node rates and node traits are not supported (see above). `nodeRates`, `nodeTraits` and `leafTrait`
   elements give a parse error with `implementation="fast"`.

## Components that required a DefaultTreeModel

Found while implementing, and changed on this branch to accept any `TreeModel`:

- `subtreeSlide` (`SubtreeSlideOperator` took a `DefaultTreeModel` to swap node rates and traits; it
  no longer does);
- `EfficientBirthDeathLikelihood` (`birthDeathLikelihood`), which only checked for a
  `DefaultTreeModel`;
- `CTMCScalePrior`'s incremental tree length (from `optimization_ctmc_prior`, merged into this branch)
  is now used with `FastTreeModel` too.

The legacy components that use node rates or node traits still require a `DefaultTreeModel`.
