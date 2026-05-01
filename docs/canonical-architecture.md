# Canonical OU Architecture

The canonical OU tree implementation is organized around a narrow, flat row-major
numerical core. Matrices in canonical hot paths should be represented as
`double[]` with index `row * dimension + column`. Use `double[][]` or EJML
objects only at legacy, BEAST, or external-library boundaries where conversion
would cost more clarity or performance than it saves.

## Module Boundaries

- `dr.evomodel.treedatalikelihood.continuous.canonical.traversal`
  owns tree state plus post-order and pre-order movement. It should not own
  gradient target math, BEAST/CDI lifecycle, or OU backend specialization.
- `dr.evomodel.treedatalikelihood.continuous.canonical.message`
  owns canonical Gaussian states, branch transition factors, local adjoints, and
  Gaussian message algebra. It should stay tree-topology agnostic.
- `dr.evomodel.treedatalikelihood.continuous.canonical.gradient`
  owns tree-level adjoint preparation and pullbacks to model gradient targets.
  It should not know BEAST/CDI lifecycle details or specialized backend classes.
- `dr.evomodel.treedatalikelihood.continuous.canonical.adapter`
  owns BEAST/CDI lifecycle, tip synchronization, root-prior adapters, dirty state,
  and cache diagnostics. It adapts external model objects into the canonical core.
- `dr.evomodel.continuous.ou.canonical`
  owns the OU kernel contracts and optional backend capability interfaces. It is
  the narrow model/kernel surface consumed by canonical tree code.
- `dr.evomodel.continuous.ou.orthogonalblockdiagonal`
  owns the private specialized orthogonal-block backend: basis caching, prepared
  branch data, transition/covariance assembly, and native gradient pullbacks.

Boundary tests in `CanonicalPackageBoundaryTest` enforce the most important
dependencies. When a new dependency violates those tests, prefer adding a small
capability interface at the boundary over importing implementation packages
directly.

## Gradient Ownership

Canonical OU tree gradients enter through
`canonical.adapter.CanonicalOUGradientAdapter`, which delegates to
`canonical.adapter.CanonicalOUTreeLikelihoodIntegrator`. The integrator owns
gradient accumulators, the cached joint-gradient values, and the reusable
`canonical.gradient.BranchGradientInputs` prepared state.

`BranchGradientInputs` is the single canonical tree-gradient preparation product:
it carries branch-local `CanonicalLocalTransitionAdjoints`, effective branch
lengths, prepared transition/basis handles from the transition cache, and the
root pre/post-order states. Selection, diffusion, stationary-mean, branch-length,
and root-mean gradients should consume this prepared state instead of preparing
their own local branch adjoints.

`BranchSpecificGradient` and `ContinuousTraitGradientForBranch` are legacy
branch-gradient machinery. They must not become dependencies of the canonical
tree-gradient backend. Canonical diagnostics may compare against legacy or finite
difference results, but those comparisons should live in tests or explicit debug
helpers rather than normal analytic dispatch.

Parallel canonical branch-gradient preparation and pullback use caller-owned
`BranchGradientWorkspace` instances, one per task-pool worker plus a reduction
workspace. Shared transition/basis data is owned by the transition cache; worker
loops should only write to worker-local scratch or to child-indexed staging slots.

## Orthogonal-Block Pullbacks

The orthogonal-block backend is the only implementation package allowed to know
about the specialized block structure. Canonical tree gradients interact with it
through model capability interfaces and caller-owned workspaces. The hot path is:

- `CanonicalSelectionGradientPullback` prepares tree-level branch adjoints.
- `SpecializedCanonicalSelectionGradientPullback` accumulates into reusable
  `SpecializedGradientWorkspace` buffers.
- `SpecializedCanonicalSelectionParameterization` and
  `CanonicalGradientPackingCapability` expose fill/accumulate APIs only.
- `orthogonalblockdiagonal` owns prepared branch/basis caches and any private
  EJML or block-specialized scratch.

Orthogonal-block pullbacks should be allocation-free after construction:
rotation, compressed-block, native-block, diffusion, and mean gradients all use
flat row-major `double[]` buffers supplied by the caller or by an explicit
workspace. Do not introduce `double[][]` into canonical gradient internals, and
do not add return-array pullback APIs to the canonical tree-gradient path.
