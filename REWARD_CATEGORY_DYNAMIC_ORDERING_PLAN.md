# Plan: reward-sorted categories with dynamic continuous-slot positioning

Branch: `feature/reward-category-dynamic-ordering` (off `rewards-aware-models`)

## Status: Phases 1-6 implemented and validated (commit `e2fb69945`)

Implemented leaner than originally sketched in ~4: rather than building a
full per-branch `EmbeddedOrdinalParameter`/cuts array, the axis itself (cut
values `0,1,...,K+1`) turns out to be branch-independent -- only the
bucket-to-category *translation* needs a per-branch insertion rank. So the
shared static-cuts decoding logic is reused unchanged, and each branch only
carries one `int` (its rank). Decode stays O(log K) per call, matching the
shared decoder's cost, not O(branchCount); refresh is O(branchCount log K),
paid once per operator call as required.

Two real, non-obvious things fell out of actually building and running this,
both worth knowing before writing XML against it:

1. **Exact-perturbation diagnostics don't work with the dynamic decoder yet.**
   `BeagleRewardDependentCtmcEdgeEvidenceProvider`'s `CategoricalRewardStateAdapter`
   dispatches on branch-rate-model *type*, not just the shared
   `RewardMixtureBranchRateModel` interface, for its diagnostic
   candidate-setting/exact-reevaluation path -- and that path is branch-index-
   independent, incompatible with a branch-specific embedding. Rather than
   silently computing wrong exact-comparison numbers, it now throws a clear
   `UnsupportedOperationException` if `dependentCtmcCompareExact` /
   `dependentCtmcDiagnostics` is used with `rewardsAwareCategoricalMixtureBranchRatesDynamic`.
   Real production sampling (the actual point of this feature) is unaffected;
   only the optional exact-comparison diagnostic path is blocked.
2. **The XML convention for "start every branch in the continuous category"
   changes.** With the static decoder, initializing `reward.category` to
   `0.5` for every branch reliably starts everything continuous, because
   category 0 is always the first axis bucket. With the dynamic decoder,
   which bucket is continuous depends on that branch's insertion rank among
   the sorted atomic rewards -- so the initial value needs to be
   `rank(initial_cts) + 0.5`, not a fixed `0.5`. Discovered by hitting an
   `-Inf` initial likelihood in the first real smoke run: with sorted atomic
   rewards `[0.0, 0.25, 0.75, 1.0]` and `cts=0.5`, rank is `2` for every
   branch, so `reward.category` needed `2.5`, not `0.5` -- see
   `tests/TestXML/testRewardsAwareRabiesHostDnaMixedDynamicSmoke.xml` for a
   worked example. Worth a clearer error message or a convenience initializer
   before this ships beyond a smoke.

Not yet done, deferred as genuine follow-up work (not blocking):
- The reward-sorted static ordering (§2) still needs to actually be applied
  in `generate_full_scale_xml.py` -- this plan only implemented the code
  side, not the XML-generation change.
- A full-scale (372-taxon) rerun comparing atomic-state reachability against
  the static-ordering baseline from `project_log.md` 2026-08-22 (only atomic
  state "Ap" was ever reached in 5 states there) has not been run yet --
  that is the actual test of whether this solves the original motivating
  problem.

## 1. Motivation

Investigated in `dep-markov/private-code/notes/project_log.md` (2026-08-21/22): the
discontinuous-HMC (DHMC) walk over an embedded categorical `reward.category`
coordinate can only cross boundaries *sequentially, one at a time, in
whatever order the categories are laid out along the coordinate axis*
(`DiscontinuousCoordinateIntegrator.traceStep` calls
`provider.getNextDiscontinuity` for the *immediately next* boundary only).
Today:

- Category 0 ("continuous") is a fixed slot `[cuts[0], cuts[1]]`, always at
  one end of the axis.
- Categories `1..K` (the atomic host states) are laid out in whatever order
  `host.rewardRates.stateIndices` happens to list them — currently arbitrary
  (alphabetical-by-species-code in the rabies XML).

Consequence: only the atomic state adjacent to the continuous slot is one
boundary-crossing away; every other atomic state requires traversing all the
states between it and the continuous slot in a single HMC step. This is a
real, previously-undiscovered reachability bias — worth fixing independently
of DHMC's actual per-crossing cost (already fixed this session).

Two independent improvements were discussed:

- **(A) Reward-sorted atomic order** — lay out categories `1..K` by reward
  value instead of arbitrary order, so "adjacent in the embedding" means
  "adjacent in reward space." This is a pure XML/data choice — see §2.
- **(B) Dynamic continuous-slot positioning** — instead of a *fixed* end
  slot, place the continuous category between whichever two (reward-sorted)
  atomic states currently bracket that branch's live `total.rewards.cts`
  value, so *every* atomic state gets a genuine one-crossing shot at being
  reached whenever `cts` drifts near it. This is the real architecture work
  this plan covers — see §3 onward.

## 2. (A) Reward-sorted atomic order — no code change, do this now

`RewardsAwareCategoricalMixtureBranchRates.getRawRewardForAtomState(stateIndex)`
maps axis position `stateIndex` (`= category - 1`) through
`host.rewardRates.stateIndices` into the concatenated
`fixedValues + varyingValues` array. Sorting the atomic axis order by reward
is entirely a matter of choosing `stateIndices`'s *values* correctly when
building the XML — no Java changes needed. This should be applied to
`runs/rabies_host_dna_smoke/scripts/generate_full_scale_xml.py`
independently of the rest of this plan, and is not blocked on anything
below.

## 3. (B) Architecture: today's global embedding vs. the target per-branch one

Confirmed by grep (`RewardMixtureCategoryDecoder`, `getCategoryCutParameter`,
`EmbeddedOrdinalParameter` usages): every consumer of the categorical
embedding today shares **one** `EmbeddedOrdinalParameter`, built from **one**
`categoryCuts` parameter, for the whole model:

| File | Usage |
|---|---|
| `RewardsAwareCategoricalMixtureBranchRates` | owns the shared `categoryDecoder`; `getUntransformedBranchRate` looks up atomic state via `categoryDecoder.getAtomicState(p)` |
| `RewardsAwareBranchModel` | owns its own shared `categoryDecoder` instance (legacy-vs-categorical branch models each construct one) |
| `RewardsAwareBranchModelGradient` | reads `rewardsAwareBranchModel.getCategoryDecoder()`, calls `isAtomic(parameterIndex)` |
| `RewardMixtureCategoricalDiscontinuousPotentialProvider` | owns a `categoryDecoder`; `getCategoryForValue`/`getNextDiscontinuity`/`getPotentialDifferenceAcrossBoundary` all route through it |
| `BeagleRewardDependentCtmcEdgeEvidenceProvider.CategoricalRewardStateAdapter` | reads `branchRates.getCategoryCutParameter()` directly to compute `representativeValueForCategory` |
| `RewardMixtureCategoricalGibbsOperator` | owns its own `categoryDecoder`, but — see §3.1 — **is not actually affected by this problem** |

`total.rewards.cts` (the parameter whose value should determine the
continuous slot's position) is **per-branch** (`dimension = branchCount`).
So "the continuous slot sits at cts's current position" cannot be a single
shared cuts array — each branch needs its own. This is the core scope driver
for everything below.

### 3.1 Gibbs operator is out of scope (good news)

`RewardMixtureCategoricalGibbsOperator` resamples a branch's category
directly from its exact full conditional (branch-local weights computed
fresh, not by boundary-walking) — it already reaches every category with
correct relative probability in one step, regardless of axis order. This
confirms the earlier session's recommendation to run it alongside DHMC as a
complementary, order-independent mixing move; it needs no changes here.

### 3.2 Design invariant (from your answer above — this drives the whole design)

> The atomic-reward order (hence the per-branch cut layout) **may change
> between MCMC iterations**, but must be **fixed for the entire duration of
> one HMC operator call** (all `nSteps`, every boundary crossing within that
> one trajectory) — never recomputed mid-walk.

This is not just a design preference — it is *required* for the DHMC
integrator to stay exact (Nishimura et al.'s reversibility/energy-exactness
proof assumes a static piecewise-constant potential landscape during one
integrated trajectory; that's also *why* `Pr(accept)=1.0` is expected for
this integrator, as seen throughout this session's validation runs). It also
directly avoids the class of bug fixed earlier this session
(`getPotentialDifferenceAcrossBoundary` redundantly recomputing state that
should have been cached once per operator call, not once per boundary
check) — the new per-branch layout must be **snapshotted once at the start
of the operator call** (the same point where `refreshEmbedding()` /
`refreshLikelihoodMessages()` are already called today) and reused
unchanged for every boundary evaluation inside that call.

## 4. Design: per-branch cut layout construction

Given, for one branch, at snapshot time:
- Sorted atomic rewards `r_1 < r_2 < ... < r_K` with associated (fixed at
  snapshot time) atomic state indices `a_1, ..., a_K` (the sort itself needs
  to be recomputed at snapshot time too, from current
  `rewardRates.getValues()` — cheap, `O(K log K)`, K is small, e.g. 17).
- This branch's current `total.rewards.cts` value `c`.

Construct:
1. `m` = number of `r_i < c` (rank of `c` among the sorted atomic rewards;
   `0 <= m <= K`).
2. A `(K+2)`-length cuts array, structurally identical in size/shape to
   today's `stateCount + 2` convention, but with the unit-width "continuous"
   slot inserted at axis position `m` instead of always at position 0:
   `cuts = [0, 1, ..., m, m+1, ..., K+1]` (same integers as today — only
   the **mapping from axis position to physical category** changes, not the
   cut values themselves).
3. Axis-position-to-category mapping for this branch, this snapshot:
   - position `< m`: atomic state `a_{position+1}` (sorted rank
     `position`, i.e. the `position`-th smallest atomic reward)
   - position `== m`: continuous
   - position `> m`: atomic state `a_{position}` (sorted rank
     `position - 1`, shifted by the inserted continuous slot)

### 4.1 The exact-tie edge case (must handle explicitly — we just spent a whole session on this)

If `c` exactly equals some `r_i`, tie-break deterministically (e.g. "ties
count as `c` being just below `r_i`," i.e. use `r_i < c` not `r_i <= c`, so
the continuous slot always sits strictly on one defined side). More
importantly: **document this loudly** and add a defensive check — the whole
day's investigation was triggered by a continuous value coinciding exactly
with an atomic reward being a genuine mathematical degeneracy in the reward
gradient, not a code bug. If `c` lands within some small epsilon of an
atomic reward, that branch's local gradient may still be numerically
delicate regardless of how the tie is broken in the embedding — this design
doesn't remove that underlying degeneracy, it only decides which side of the
boundary to put the continuous slot on. Flag this as a known residual risk
in the tests (§6) rather than something this plan claims to fully solve.

## 5. Implementation phases (new, opt-in classes only — per your answer above)

Nothing below modifies `RewardMixtureCategoryDecoder`,
`RewardsAwareCategoricalMixtureBranchRates`, `RewardsAwareBranchModel`,
`RewardMixtureCategoricalDiscontinuousPotentialProvider`, or
`BeagleRewardDependentCtmcEdgeEvidenceProvider` in place. Every existing
validated XML/test keeps using the current shared-embedding classes
unchanged.

### Phase 1 — `PerBranchCategoryLayout` (new, pure logic, no BEAST model wiring)

A small, standalone, thoroughly-unit-testable class implementing exactly the
construction in §4: given sorted atomic rewards + one `cts` value, produce
the axis-position-to-category mapping for one branch. No `Parameter`/
`Model`/tree dependencies — easy to test exhaustively (rank boundaries,
`m=0`, `m=K`, ties, duplicate reward values, `K=1`).

### Phase 2 — `PerBranchRewardMixtureCategoryDecoder` (new, parallel to `RewardMixtureCategoryDecoder`)

Wraps one `PerBranchCategoryLayout` **per branch**, built from:
`rewardRates.getValues()` (for the sort) and `ctsParameter` (per-branch
values). Exposes the same surface `RewardMixtureCategoryDecoder` does
(`getCategoryForParameterIndex`, `isAtomic`, `getAtomicState`,
`getLowerCut`/`getUpperCut`, `getNextBoundary`) but routes every call
through the *branch-specific* layout for that `parameterIndex`.
`refreshEmbedding()` becomes "rebuild all `branchCount` layouts" — O(branch
count × K log K), cheap at real scale (742 × 17 log 17 is trivial;
contrast with the O(branchCount) *redundant* work bug fixed earlier this
session, which was the *same* asymptotic class of cost but paid on every
single boundary crossing instead of once per operator call — the §3.2
invariant is exactly what keeps this new cost paid once, not repeatedly).

### Phase 3 — `RewardsAwareCategoricalMixtureBranchRatesDynamic` (new branch-rate model)

Copy of `RewardsAwareCategoricalMixtureBranchRates` (§3's first row) using
`PerBranchRewardMixtureCategoryDecoder` in place of the shared decoder.
Same `RewardMixtureBranchRateModel` interface, so it's a drop-in
alternative wherever `RewardsAwareCategoricalMixtureBranchRates` is
currently accepted (`RewardsAwareBranchModel`,
`BeagleRewardDependentCtmcEdgeEvidenceProvider`'s
`RewardMixtureBranchRateModel`-typed constructor parameter, etc.) — check
this interface-level substitutability carefully before assuming zero
further changes are needed downstream; some downstream code may
`instanceof`-check for the concrete class rather than the interface (grep
before starting).

### Phase 4 — `RewardMixtureCategoricalDiscontinuousPotentialProviderDynamic` (new, or generalize via the interface)

The existing `RewardMixtureCategoricalDiscontinuousPotentialProvider`
constructs its own internal `categoryDecoder`
(`RewardMixtureCategoryDecoder`) directly rather than taking one as a
dependency — needs either a parallel `...Dynamic` class, or (preferred if
feasible) refactor the existing class to depend on a small interface both
`RewardMixtureCategoryDecoder` and `PerBranchRewardMixtureCategoryDecoder`
implement, then construct the right concrete decoder at the parser level.
Prefer the interface-extraction approach if it doesn't require touching
existing call sites' behavior — reduces duplicate boundary-crossing/weight
logic between the two variants, which is exactly the kind of logic that's
already had one serious correctness bug this session (the
`getPotentialDifferenceAcrossBoundary` fix) and shouldn't be forked
carelessly.

### Phase 5 — XML wiring

New parser(s), or a boolean attribute
(`dynamicContinuousPosition="true"`) recognized by a small dedicated parser
class rather than overloading the existing, validated
`RewardsAwareCategoricalMixtureBranchRatesParser` — matches this project's
established pattern of adding new dev-only parsers alongside existing ones
(e.g. how the categorical DHMC layer itself was added in June without
touching the legacy indicator/atom parsers).

### Phase 6 — Testing

- Unit tests for `PerBranchCategoryLayout` (Phase 1) — pure, exhaustive,
  no BEAST scaffolding needed; this is where most of the correctness
  confidence should come from, since it's the one piece with no BEAGLE/tree
  dependencies.
- A regression test mirroring
  `RewardMixtureCategoricalDiscontinuousPotentialProviderTest`'s existing
  manual/exact/BEAGLE-delta-agreement pattern, but with the dynamic decoder
  — confirm evidence/weights still match exact perturbation to roundoff.
- A small XML smoke analogous to
  `tests/TestXML/testRewardsAwareRabiesHostDnaMixedSmoke.xml` (5 real
  taxa), using the new dynamic branch-rate model, checked with a real
  `gradientCheckCount` pass (informed by today's finding: initialize `cts`
  off any exact atomic-reward coincidence, e.g. `0.5001`-style offset, not
  `0.5`).
- Full-scale rerun (`runs/rabies_host_dna_smoke/scripts/generate_full_scale_xml.py`,
  extended to emit the dynamic wiring) as the actual test of the original
  motivating question: does a longer chain now visit atomic states that
  were unreachable under the static ordering baseline recorded in
  `project_log.md` 2026-08-22 (only atomic state 1, "Ap," was ever reached
  in 5 states under static ordering)?

## 6. Explicit non-goals for this branch (V1)

- Handling reordering *within* a single HMC trajectory — not just
  deferred, actually prohibited by the correctness invariant in §3.2.
- Automatically detecting/relabeling permutation ambiguity if atomic reward
  values cross (that's the `OneZeroOneShuffleGibbsOperator` problem in a
  continuous-reward-value setting — a related but separate feature, not
  bundled here).
- Migrating any existing validated XML/run to the new classes — this stays
  strictly additive until someone deliberately opts a new run into it.

## 7. Suggested order of work

1. Phase 1 (`PerBranchCategoryLayout` + its unit tests) — no BEAST
   dependencies, fastest to get right and build confidence in the core
   algorithm including the tie-handling edge case.
2. Phase 2 (`PerBranchRewardMixtureCategoryDecoder`) with its own focused
   unit tests against a hand-built tiny `Parameter` fixture (no tree/
   likelihood needed yet).
3. Phases 3-4 together (the two new model/potential classes), validated via
   the small 5-taxon XML smoke (Phase 6) before touching full scale.
4. Phase 5 (XML parser) once 3-4 are validated in Java-driver tests, so the
   parser wiring is the last, lowest-risk step.
5. Full-scale rerun (Phase 6 tail) as the final validation, comparing
   atomic-state reachability against today's static-order baseline.
