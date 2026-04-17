# CLAUDE.md — beast-mcmc-time-series-e0

## Project overview

This branch is rewriting the continuous-trait tree likelihood around a **canonical-form OU pathway**.

The current architectural center is no longer the older moment-form tree passer from `CLAUDEOld.md`. The main production story is now:

- **OU first**
- **canonical form first**
- **tree-native message passing**
- **exact gradients for HMC**

The canonical pathway is intended to be the new reference architecture for OU on trees. The older moment-form tree framework still exists in the repository, but it is no longer the conceptual source of truth for the new OU design.

---

## Build and test

```bash
ant compile-all
ant junit-canonical
```

Useful focused test commands:

```bash
java -cp "build:lib/*" junit.textui.TestRunner \
  test.dr.evomodel.treedatalikelihood.continuous.CanonicalOrthogonalBlockFrechetExactTest

java -cp "build:lib/*" junit.textui.TestRunner \
  test.dr.evomodel.treedatalikelihood.continuous.CanonicalGaussianMessagePasserOUGradientTest

java -cp "build:lib/*" junit.textui.TestRunner \
  test.dr.evomodel.treedatalikelihood.continuous.CanonicalOUXmlWiringTest
```

## Running XMLs

When running BEAST XMLs from this branch, always disable citation-file generation:

```bash
java -cp "build:lib/*" dr.app.beast.BeastMain -citations_off <other args> <file.xml>
```

Keep `-citations_off` on by default for smoke runs, debugging runs, static checks, and real-data XML checks.

---

## Canonical architecture

The canonical OU path is split across two layers:

1. **Tree-side canonical message passing**
2. **Time-series canonical Gaussian algebra**

The tree side owns traversal, root/tip wiring, and branch-local gradient orchestration. The time-series side owns canonical Gaussian states, transitions, and low-level canonical operations.

### Canonical tree contracts

All canonical tree interfaces live in:

`src/dr/evomodel/treedatalikelihood/continuous/framework/`

Core classes:

- `CanonicalTreeMessagePasser`
  - canonical post-order / pre-order / gradients / store-restore API
- `CanonicalBranchTransitionProvider`
  - maps child node index to an exact canonical branch factor
- `CanonicalRootPrior`
  - canonical root seeding and canonical root likelihood integration
- `CanonicalTipObservation`
  - exact observed/missing tip representation for the canonical path

Important design decision:

- canonical tree messages are stored directly as `CanonicalGaussianState`
- canonical branch factors are stored directly as `CanonicalGaussianTransition`
- the canonical OU path does **not** wrap these in a second tree-specific math type

### Tree-side implementation

The main tree implementation is:

- `SequentialCanonicalOUMessagePasser`

This is the reference canonical OU tree algorithm. It currently handles:

- post-order canonical likelihood propagation
- pre-order canonical propagation
- exact observed tips
- partially missing tips
- `Q`, `A`, `mu`, and branch-length gradients
- MCMC `storeState` / `restoreState` / `acceptState`

The main end-to-end wiring class is:

- `CanonicalOUMessagePasserComputer`

This is the canonical analogue of the older message-passer computer classes. It owns:

- the canonical passer
- the canonical branch-transition provider
- the canonical root prior
- loading canonical tip observations

### Adapters and providers

The current canonical OU adapters live in:

`src/dr/evomodel/treedatalikelihood/continuous/adapter/`

Most important:

- `HomogeneousCanonicalOUBranchTransitionProvider`
  - homogeneous OU branch provider
  - delegates branch construction to `OUProcessModel.fillCanonicalTransition(...)`
  - snapshots diffusion precision -> covariance lazily
  - caches canonical branch transitions per node
- `CanonicalConjugateRootPriorAdapter`
  - canonical adapter over `ConjugateRootTraitPrior`
  - integrates the root message in canonical form
- `CanonicalTipObservationAdapter`
  - XML/data-model boundary adapter into `CanonicalTipObservation`

Current gradient support is intentionally restricted:

- gradients currently require `HomogeneousCanonicalOUBranchTransitionProvider`
- heterogeneous canonical providers throw a clear `Not yet implemented` error

---

## Message-passing algorithm

### Stored states

The canonical passer stores three state families:

- `postOrder[node]`
  - upward canonical message at each node
- `preOrder[node]`
  - downward canonical message at each node
- `branchAboveParent[child]`
  - the explicit message just above a branch, before pushing through that branch

This last one is important:

- the production path uses the **explicit `branchAboveParent` message**
- it is stored directly during pre-order traversal
- it is the message used later in branch-local gradient assembly
- recovery-from-child-above is not the production architecture

### Post-order

For each internal node:

1. build each child-to-parent canonical message
2. combine child messages canonically
3. store the result in `postOrder[parent]`

For tips:

- fully observed tips are handled as exact observation factors
- partially missing tips use missing-safe precision/variance algebra, but remain in the canonical tree pathway

At the root:

- the upward root message is integrated against `CanonicalRootPrior`

### Pre-order

At the root:

- if the root is fixed, the passer stores the fixed root value explicitly
- otherwise, the root prior seeds a canonical root message

For each child:

1. combine the parent's downward message with sibling information
2. store that explicit parent-side message in `branchAboveParent[child]`
3. push it through the child branch in canonical form
4. store the result in `preOrder[child]`

This is the downward production path used by the gradient code.

---

## Canonical Gaussian backend

The canonical tree path relies on the time-series canonical representation layer in:

`src/dr/inference/timeseries/representation/`

Key types:

- `CanonicalGaussianState`
- `CanonicalGaussianTransition`
- `CanonicalGaussianUtils`

And on the canonical Gaussian operations in:

- `CanonicalGaussianMessageOps`

This utility now owns the shared canonical algebra used by both the tree passer and the canonical Kalman engines:

- combining canonical states
- forward / backward push through canonical transitions
- pair-posterior assembly
- canonical normalization shifts
- block marginalization utilities

The canonical path no longer depends architecturally on `KalmanLikelihoodEngine` as its math host. Shared dense linear algebra was extracted into:

- `GaussianMatrixOps`

---

## OU parametrization and adjoint defaults

For OU, the intended default chart is the **orthogonal block-diagonal chart**.

That means:

- the default OU XML path should be interpreted as orthogonal-block unless dense is explicitly requested
- the preferred adjoint machinery is the **stationary Lyapunov / orthogonal-block path**
- the dense chart remains available for generality, but it is not the intended default production route

Parser policy now reflects that:

- orthogonal block is the implicit default
- dense charts must be requested explicitly with `selectionChart="dense"`

Relevant classes:

- `OrthogonalBlockDiagonalSelectionMatrixParameterization`
- `OUProcessModel`
- `OUSelectionChartParserHelper`

---

## Exact orthogonal-block Fréchet backend

One important recent architectural change is the removal of numerical quadrature and generic tiny matrix exponentials from the hot orthogonal-block Fréchet production path.

The exact block Fréchet plan lives in:

- `BlockDiagonalFrechetExactPlan`

It now uses:

- closed-form `1x1`, `1x2`, and `2x1` block formulas
- a specialized **equal-diagonal `2x2-2x2`** path
- analytic small-root series stabilization near the nilpotent regime
- generic exact fallback only when a `2x2` block pair is not equal-diagonal

Why this matters:

- the orthogonal polar chart produces `2x2` blocks of the form
  - `[a  u; l  a]`
- so each block has equal diagonals
- that lets the Fréchet map be reduced to the basis
  - `{X, N_L X, X N_R, N_L X N_R}`
- which is much cheaper than solving the generic `4x4` Sylvester problem for every block pair

This optimization is now live in the Java path and validated against dense fallback tests.

Supporting backprop classes:

- `BlockDiagonalFrechetHelper`
- `BlockDiagonalExpSolver`
- `BlockDiagonalLyapunovSolver`
- `BlockDiagonalLyapunovAdjointHelper`

---

## Gradient architecture

The canonical gradient story is:

1. run post-order
2. run pre-order
3. for each branch, build a canonical branch-local contribution from:
   - `postOrder[child]`
   - `branchAboveParent[child]`
   - the canonical branch transition
4. pull local adjoints back to model parameters

Main branch-local canonical adjoint types:

- `CanonicalBranchMessageContribution`
- `CanonicalLocalTransitionAdjoints`
- `CanonicalTransitionAdjointUtils`

### `Q`

`computeGradientQ(...)` works canonically end to end, including the conjugate-root diffusion-scaling term from `CanonicalRootPrior.getDiffusionScale()`.

For the orthogonal-block chart, the diffusion pullback uses the exact native path rather than forcing everything through a dense fallback.

### `A`

`computeGradientA(...)` supports:

- generic dense selection gradients
- orthogonal-block native gradients

For orthogonal block, the production path is the native pullback through:

- `OrthogonalBlockDiagonalSelectionMatrixParameterization`

### `mu`

`computeGradientMu(...)` is implemented canonically and includes the canonical root contribution when appropriate.

### branch lengths

`computeGradientBranchLengths(...)` is implemented for the canonical OU path.

At the moment this uses a local branch finite-difference route, not a final closed-form branch-time pullback.

---

## XML wiring

The canonical XML entry point is:

```xml
<traitDataLikelihood implementation="canonical" ... />
```

Current parser behavior:

- canonical `traitDataLikelihood` currently supports **OU only**
- canonical `traitDataLikelihood` does **not** support `reconstructTraits="true"`

Canonical gradient XML is wired for:

- `meanGradient implementation="canonical"`
- `precisionGradient implementation="canonical"`
- `attenuationGradient implementation="canonical"`

These canonical gradient parsers require:

- `traitDataLikelihood implementation="canonical"`

The canonical XML path has been validated on the orthogonal-block OU real-data XMLs already used in this branch.

---

## Tests

Canonical tests are intentionally split from the older moment-form passer tests.

Canonical-focused tests include:

- `CanonicalGaussianMessagePasserOULikelihoodTest`
- `CanonicalGaussianMessagePasserOUGradientTest`
- `CanonicalOrthogonalBlockFrechetExactTest`
- `CanonicalOUXmlWiringTest`

Other useful validation tests:

- `OUDiffusionKernelBridgeValidationTest`
- `OUProcessModelTest`

Run them together with:

```bash
ant junit-canonical
```

---

## Current limitations

- canonical tree likelihood currently targets **OU**, not BM
- canonical XML currently requires `traitDataLikelihood implementation="canonical"` and OU
- `reconstructTraits="true"` is not supported in canonical mode
- gradients currently support only `HomogeneousCanonicalOUBranchTransitionProvider`
- exact observed/missing tip semantics are supported; finite-noise tip observations are not yet the canonical API
- heterogeneous canonical gradient providers are intentionally `Not yet implemented`

---

## Practical guidance

If you are working on the new OU path, start from these files:

- `SequentialCanonicalOUMessagePasser`
- `CanonicalOUMessagePasserComputer`
- `HomogeneousCanonicalOUBranchTransitionProvider`
- `CanonicalConjugateRootPriorAdapter`
- `CanonicalGaussianMessageOps`
- `OrthogonalBlockDiagonalSelectionMatrixParameterization`
- `BlockDiagonalFrechetExactPlan`

If you are debugging correctness:

1. check canonical XML wiring first
2. check `branchAboveParent` storage and usage
3. compare orthogonal native gradients against dense fallback
4. compare the exact Fréchet path against `CanonicalOrthogonalBlockFrechetExactTest`

If you are debugging performance:

1. profile branch-transition caching in `HomogeneousCanonicalOUBranchTransitionProvider`
2. profile canonical Gaussian linear algebra in `CanonicalGaussianMessageOps`
3. profile orthogonal-block Fréchet / Lyapunov helpers
4. prefer removing recomputation over micro-optimizing scalar transcendentals

