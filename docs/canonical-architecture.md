# Canonical OU Architecture

The canonical OU tree implementation is organized around a narrow, flat row-major
numerical core. Matrices in canonical hot paths should be represented as
`double[]` with index `row * dimension + column`. Use `double[][]` or EJML
objects only at legacy, BEAST, or external-library boundaries where conversion
would cost more clarity or performance than it saves.

## Delegate Split: Legacy vs. `continuousNew`

Two `DataLikelihoodDelegate` implementations exist side by side, both behind the
shared `dr.evomodel.treedatalikelihood.continuousNew.CtsCompatibility`
interface. Generic code (`BranchRateGradient`, `BranchConditionalDistributionDelegate`,
`CanonicalSelectionParameterGradient`, `CanonicalMeanParameterGradient`,
`CanonicalPrecisionGradient`, `BranchRateGradientParser`, `NodeHeightGradientParser`,
`ConditionalTraitSimulationHelper`, etc.) should depend on `CtsCompatibility`,
not on either concrete class, unless it specifically needs legacy-only or
new-only behavior.

- `dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate`
  (legacy) only ever runs the legacy CDI/SafeIntegration path
  (`calculateLegacyCdiLikelihood`). It never touches `likelihoodEngine`. Its
  `enable*Likelihood()` methods (`enableCanonicalLikelihood`,
  `enableCdiBMLikelihood`, `enableMomentBMLikelihood`,
  `enableCanonicalOULikelihood`) all throw `UnsupportedOperationException` by
  design — this class cannot run any of the newer backends, on purpose.
- `dr.evomodel.treedatalikelihood.continuousNew.ContinuousDataLikelihoodDelegateNew`
  is the eventual replacement. It never builds a real legacy
  `ContinuousDiffusionIntegrator`: `cdi` stays `null` unless constructed with
  `initializeLegacyCdi=true`, and that branch of its constructor is itself
  deprecated (it unconditionally throws). All computation instead flows
  through `likelihoodEngine` (a `ContinuousTreeLikelihoodEngine`), selected via
  `ContinuousTreeLikelihoodEngineFactory`: `momentBM` (native moment-form BM,
  the current BM production backend), `canonical` (dispatches to `canonicalOU`
  for OU models), `canonicalOU`. `cdiBM`/`CdiBMContinuousTreeLikelihoodEngine`
  and the plain legacy-CDI engine selection are commented-out/dead — the CDI-
  backed BM engine has been superseded by the native moment-form engine, so
  `enableCdiBMLikelihood()` is currently a no-op and any implementation mode
  that calls it silently keeps whatever engine the constructor already
  selected.

`dr.evomodelxml.treedatalikelihood.ContinuousDataLikelihoodParser` is the only
production construction site for either class: `implementation="legacy"`
builds the old delegate (with `initializeLegacyCdi` computed from the
implementation mode); every other `implementation=` value builds
`ContinuousDataLikelihoodDelegateNew` (always with `initializeLegacyCdi=false`)
and then calls `applyImplementationMode(...)`, which invokes the matching
`enable*Likelihood()`. The constructor's own baseline engine selection is
gated to homogeneous-BM diffusion models only — for OU/drift models it leaves
`likelihoodEngine` unset and relies entirely on the `enable*Likelihood()` call
that always follows immediately.

**Known limitation:** `ContinuousDataLikelihoodDelegateNew.getReport()` calls
`MultivariateTraitDebugUtilities.getTreeDrift`/`getGraphDrift`, which for a
genuine drift/OU diffusion model needs a real `cdi` to read cached branch
displacements — something the new delegate never builds. Calling `getReport()`
(e.g. via a `<report>` XML element) on an OU/drift model through the new
delegate throws `NullPointerException`. Homogeneous BM avoids this because
`HomogeneousDiffusionModelDelegate.getAccumulativeDrift` ignores `cdi`
entirely and just returns the prior mean.

## Module Boundaries

**Still under `dr.evomodel.treedatalikelihood.continuousNew`** (tree-dependent,
or deliberately kept here — see "Tree-vs-Math Separation" below):

- `continuousNew` (root) holds the delegate (`ContinuousDataLikelihoodDelegateNew`),
  the `CtsCompatibility` contract, the engine factory/interface, and tree-side
  branch-process/transition-provider glue (`ContinuousBranchProcessModel`,
  `HomogeneousContinuousBranchProcessModel`,
  `ContinuousBranchTimeModel`/`DefaultContinuousBranchTimeModel`,
  `GaussianBranchTransitionProvider`/`KernelBasedGaussianBranchTransitionProvider`/
  `OUGaussianBranchTransitionProvider`).
- `continuousNew.message` holds the small post-order/pre-order message-provider
  marker interfaces (`ContinuousMessageProvider`,
  `ContinuousCanonicalPostOrderMessageProvider`,
  `ContinuousCanonicalPreOrderMessageProvider`) — distinct from
  `canonical.message` below, which owns the actual Gaussian message algebra.
- `continuousNew.observationmodel` owns tip observation models and their
  workspaces (`CanonicalTipObservationModel`, `IdentityCanonicalTipObservationModel`,
  `PartialIdentityTipProjection`, `TipObservationPartition`, etc.) plus the
  narrow `CanonicalTipBranchMessageProvider` interface that
  `CanonicalNativeBranchMessageProvider` extends. Not tree-traversal-dependent
  itself, but kept here deliberately rather than relocated (see below).
- `continuousNew.canonical.traversal` owns tree state plus post-order and
  pre-order movement. It should not own gradient target math, BEAST/CDI
  lifecycle, or OU backend specialization.
- `continuousNew.canonical.message` owns canonical Gaussian states, branch
  transition factors, local adjoints, and Gaussian message algebra. Tree-topology
  agnostic in principle (and shared with the time-series Kalman engine — see
  below) but kept here deliberately rather than relocated.
- `continuousNew.canonical.gradient` owns tree-level adjoint preparation and
  pullbacks to model gradient targets. It should not know BEAST/CDI lifecycle
  details or specialized backend classes. Most of this package is pure
  adjoint/pullback algebra with no tree dependency, but a small number of files
  (the engines that walk the tree to assemble `BranchGradientInputs`,
  `StateStoreBranchGradientMessageProvider`) are genuinely tree-dependent;
  splitting the two apart is tracked but not yet done (see below).
- `continuousNew.canonical.adapter` owns BEAST/CDI lifecycle, tip
  synchronization, root-prior adapters, dirty state, and cache diagnostics. It
  adapts external model objects into the canonical core.

**Relocated under `dr.evomodel.continuous`** (pure math, no tree dependency;
also used outside the tree-likelihood subsystem — see below):

- `dr.evomodel.continuous.canonical` (root) holds `BMProcessModel` and
  `ContinuousProcessModel`: branch-local BM process math that fills
  caller-owned buffers given an effective branch time.
- `dr.evomodel.continuous.canonical.math` owns pure flat-array linear algebra
  and Gaussian-conversion utilities (`MatrixOps`, `GaussianFormConverter`,
  `CanonicalInversionResult`). New canonical hot-path matrix code should use
  this package directly instead of `double[][]`/EJML.
- `dr.evomodel.continuous.canonical.contribution` owns branch-local canonical
  contribution strategies (exact tip, partial tip, Gaussian-link tip, internal
  node, missing, fixed-root) assembled from a `BranchGradientMessageBundle`.
- `dr.evomodel.continuous.canonical.workspace` and `.scheduling` own reusable
  buffers/workspace factories and parallel chunk-sizing policies, respectively.
- `dr.evomodel.continuous.ou.canonical.branch` owns OU-specific branch algebra
  (`OUCanonicalTransitionBuilder`, `OUCanonicalBranchGradientBridge`, etc.) —
  nested under the pre-existing `ou.canonical` (see below) rather than the
  process-agnostic `dr.evomodel.continuous.canonical` parent, since this one
  genuinely is OU-specific.

**Pre-existing, unaffected by the split above:**

- `dr.evomodel.continuous.ou.canonical` owns the OU kernel contracts and
  optional backend capability interfaces. It is the narrow model/kernel
  surface consumed by canonical tree code.
- `dr.evomodel.continuous.ou.diagonal` owns the diagonal-selection OU
  kernel/parameterization used by the native diagonal fast path
  (`DiagonalOUBranchKernel`, `DiagonalSelectionMatrixParameterization`).
- `dr.evomodel.continuous.ou.blockdiagonal` owns the private specialized
  block-diagonal backend: basis caching, prepared branch data,
  transition/covariance assembly, and native gradient pullbacks, plus its own
  `math` subpackage for the block-space Fréchet/Lyapunov solvers.

Boundary tests in `CanonicalPackageBoundaryTest` enforce the most important
dependencies (e.g. the canonical tree layer must not import the block-diagonal
backend, canonical machinery must not use wildcard imports/direct
`System.err`/`System.out`/`double[][]`, block-diagonal buffers must be
allocated only at construction time, and the relocated `dr.evomodel.continuous.canonical`
packages must not import back into `canonical.traversal`/`canonical.adapter`
or `dr.evolution.tree`). When a new dependency violates those tests, prefer
adding a small capability interface at the boundary over importing
implementation packages directly.

## Tree-vs-Math Separation

`private-code/docs/continuousNew_tree_math_separation_plan.md` tracks this
split in full, including the consumer-dependency evidence behind each
package's target location. Status:

- **Moved** (Phase 2): `canonical.math`, `canonical.contribution`,
  `canonical.workspace`, `canonical.scheduling`, `BMProcessModel`/
  `ContinuousProcessModel`, and `canonical.branch` (the one OU-specific
  package, nested under `ou.canonical.branch` rather than the shared parent)
  — see "Module Boundaries" above for their new homes.
- **Deliberately kept in `continuousNew`**: `canonical.message` and
  `observationmodel`. Both are pure math with no tree dependency (and
  `canonical.message` is additionally shared with the
  `dr.inference.timeseries` Kalman engine), so they were candidates for the
  same move, but the decision was made to leave them in place for now.
- **Not yet done** (Phase 3): `canonical.gradient` is mixed — most of it is
  pure adjoint/pullback algebra (also shared with `dr.inference.timeseries`),
  but a few files genuinely walk the tree. Splitting it file-by-file is
  tracked but not started.

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

## Block-Diagonal Pullbacks

The block-diagonal backend is the only implementation package allowed to know
about the specialized block structure. Orthogonal charts live there as a
special case of the general `A = R D R^{-1}` pathway. Canonical tree gradients
interact with the backend through model capability interfaces and caller-owned
workspaces. The hot path is:

- `CanonicalSelectionGradientPullback` prepares tree-level branch adjoints.
- `SpecializedCanonicalSelectionGradientPullback` accumulates into reusable
  `SpecializedGradientWorkspace` buffers.
- `SpecializedCanonicalSelectionParameterization` and
  `CanonicalGradientPackingCapability` expose fill/accumulate APIs only.
- `blockdiagonal` owns prepared branch/basis caches and any private
  EJML or block-specialized scratch.

Block-diagonal pullbacks should be allocation-free after construction:
rotation, compressed-block, native-block, diffusion, and mean gradients all use
flat row-major `double[]` buffers supplied by the caller or by an explicit
workspace. Do not introduce `double[][]` into canonical gradient internals, and
do not add return-array pullback APIs to the canonical tree-gradient path.
