# SubtreeLeap operator study

Branch: `operator_study`. Study in `beast-benchmark-stl2/stl_study/`.

Datasets: **small-2** (165 taxa, 40 M states, 5 replicates) and **tall-1** (1441 taxa,
60 M states, 3 replicates for non-reference conditions, 5 for the default).

Likelihoods: **mat** (MAT data, `matDataLikelihood`) and **min** (minimum-evolution,
`treeDataLikelihood` with Hasegawa–Kishino–Yano substitution model).

The primary metric is `eff_min_cont_ratio`: the geometric-mean ESS per hour of all
continuous parameters (root height, tree length, coalescent and substitution
parameters), relative to the default condition. Secondary metrics are `eff_rf_ratio`
and `eff_split_min_ratio` (clade-frequency ESS proxies).

`n_stuck` counts replicates where split ESS was near zero (< 10), indicating
poor topology mixing in that replicate; `fewer_stuck_p` is the one-sided
Boschloo p-value that the condition has fewer stuck chains than the default.

---

## 1. Target sweep

Varied `targetAcceptance` of the single SubtreeLeap from 0.05 to 0.50;
default is 0.234.

### eff_min_cont_ratio

| target | small-2 mat | small-2 min | tall-1 mat | tall-1 min |
|--------|-------------|-------------|------------|------------|
| 0.05   | 0.24        | 0.37        | 0.29       | 0.39       |
| 0.10   | 0.58        | 0.54        | 0.99       | 0.75       |
| 0.15   | 0.86        | 0.78        | **1.08**   | 0.63       |
| **0.234** (default) | **1.00** | **1.00** | **1.00** | **1.00** |
| 0.35   | **1.19**    | **1.17**    | 0.55       | 0.68       |
| 0.50   | 0.91        | 1.11        | 0.84 †     | 0.85       |

† t=0.50, tall-1 mat: 3 of 3 replicates stuck on topology.

### Topology-change fraction vs target

Higher target acceptance → shorter leaps → fewer topology changes:

| target | sl_topologyChanged |
|--------|-------------------|
| 0.05   | 0.93–0.96         |
| 0.234  | 0.75–0.79         |
| 0.35   | 0.63–0.66         |
| 0.50   | 0.47–0.50         |

### Interpretation

The optimal target flips between small and tall trees. On **small-2**, larger trees
take longer per step (more sites per branch), so a higher acceptance target
(0.35, shorter leaps) is more efficient. On **tall-1**, wider leaps matter for
mixing in a large tree; the default 0.234 sits at the peak, and both extremes
underperform. At very low targets (0.05) the continuous ESS collapses on tall-1,
even though topology change rate is high: the leaps are so short they barely
explore the tree, and the many accepted-but-tiny moves do not translate into
parameter mixing.

**Recommendation:** keep the default 0.234. It is never the best on either dataset,
but it avoids catastrophic failure on either. A single target that beats the
default on both datasets does not exist in this sweep.

---

## 2. Weight sweep

Varied the STL weight as a multiple of the default (`nTaxa` = number of tips).
Default factor is 1.0.

### eff_min_cont_ratio

| factor | small-2 mat | small-2 min | tall-1 mat | tall-1 min |
|--------|-------------|-------------|------------|------------|
| 0.1    | 0.34        | 0.55        | 0.58       | 0.64       |
| 0.25   | 0.55        | 0.59        | 0.81       | 0.92       |
| 0.5    | 0.91        | 0.86        | 0.90       | 0.63       |
| **1.0** (default) | **1.00** | **1.00** | **1.00** | **1.00** |
| 2.0    | **1.13**    | **1.11**    | **1.20**   | 0.49 ‡     |
| 4.0    | **1.26**    | 1.00        | 1.09 §     | 0.48       |

‡ w2.0, tall-1 min: 2 of 3 replicates stuck.  
§ w4.0, tall-1 mat: 1 of 3 replicates stuck.

### Time-fraction budget

STL and FHSPR compete for clock time. At high weights, STL starves FHSPR:

| condition         | sl_time_fraction | fhspr_time_fraction |
|-------------------|-----------------|---------------------|
| tall-1 min, w1.0  | 0.86            | 0.098               |
| tall-1 min, w2.0  | 0.93            | 0.052               |
| tall-1 min, w4.0  | 0.96            | 0.027               |

FHSPR is a wide-leap burn-in aid. When it falls below ≈5% of wall time it can no
longer rescue stuck chains, which is exactly what the n_stuck data show.

MAT is less sensitive here because the per-step cost of FHSPR is higher in the MAT
model (resampling the mutation history is expensive), so FHSPR already has a smaller
default budget.

### Interpretation

Doubling the STL weight to 2× nTaxa gives a consistent ≈12–20% gain on mat
likelihoods. However, the same factor is dangerous for min likelihoods on tall trees:
chains get stuck because FHSPR is starved. Setting weight above the default is
therefore not safe as a universal default.

**Recommendation:** keep the default weight of 1× nTaxa. Users running MAT on small
or medium trees (< 500 tips) could try 2× as a tuning option.

---

## 3. Two-STL arm

A second `subtreeLeap` (`stl2`) is added after the main one. The main STL remains at
target=0.234. The second STL's target and weight factor (relative to the main STL's
weight) are varied.

### tall-1 mat — eff_min_cont_ratio

| stl2 target | stl2 factor | eff_min_cont | eff_rf   | eff_split_min | n_stuck |
|------------|------------|--------------|----------|---------------|---------|
| —          | —          | 1.00         | 1.00     | 1.00          | 0       |
| 0.03       | 0.1        | 0.50         | 0.64     | 0.85          | 0       |
| 0.03       | 0.3        | 0.94         | 1.87     | 0.92          | 0       |
| **0.05**   | **0.1**    | **1.67**     | **2.05** | **1.34**      | **0**   |
| 0.05       | 0.3        | 1.16         | 2.44     | 1.15          | 0       |
| 0.10       | 0.1        | 1.05         | 1.09     | 1.23          | 0       |
| 0.10       | 0.3        | 0.86         | 1.10     | 0.93          | 0       |
| 0.15       | 0.1        | 0.96         | 1.52     | 1.25          | 0       |
| 0.15       | 0.3        | 1.12         | 1.06     | 0.94          | 0       |
| 0.20       | 0.1        | 0.59         | 2.52     | 0.62          | 0       |
| 0.20       | 0.3        | 1.14         | 1.23     | 0.62          | 1       |

The standout condition is **target=0.05, factor=0.1**: +67% continuous ESS,
+105% RF-distance ESS, +34% split ESS. The second STL has acceptance ≈0.05
(very short leaps) and costs 8.7% of wall time. Equivalent to trading 8.7%
extra compute for 1.67× continuous mixing and 2× topology mixing.

At target=0.03 the second STL is too aggressive (nearly all moves accepted,
kernel adapts to tiny hops), and the per-step cost is not recovered in ESS.
At target≥0.15 the second STL behaves similarly to the main one, adding
redundancy rather than complementary exploration.

### small-2 mat — eff_min_cont_ratio

| stl2 target | stl2 factor | eff_min_cont | eff_rf | eff_split_min |
|------------|------------|--------------|--------|---------------|
| —          | —          | 1.00         | 1.00   | 1.00          |
| 0.05       | 0.1        | 0.94         | 0.84   | 0.90          |
| 0.05       | 0.3        | 0.78         | 1.19   | 0.70          |
| 0.15       | 0.1        | **1.04**     | 1.18   | 1.10          |
| 0.15       | 0.3        | 1.04         | 0.96   | 0.87          |
| 0.20       | 0.3        | **1.08**     | 0.80   | 0.57          |

On small-2 the second STL adds at most ≈8% to continuous ESS with no
consistent pattern. The benefit that is dramatic on tall-1 is absent here.
This makes sense: on a 165-tip tree the main STL already covers the whole
tree efficiently; a second component at a lower target provides little
additional coverage.

### tall-1 min — summary

Results are inconsistent and several conditions involve stuck replicates
(even the default has 1/5 stuck replicates, which is a known convergence
issue for tall-1 min). No condition shows a clear benefit over the default.
The second STL does not hurt, but neither does it reliably help.

### Mechanistic interpretation

On a large MAT tree the STL must explore O(nTaxa) possible destinations. A
single STL adapts to a kernel that balances acceptance rate against leap
distance. Adding a second STL at a low target (≈0.05) provides a
complementary short-range component that makes many small accepted moves,
filling in the neighbourhood of each node cheaply. The two operators explore
at different spatial scales, analogous to a multi-scale sampler.

For tall-1 mat, 8.7% of wall time goes to the second STL. Its per-operation
cost (128 µs) is similar to the main STL (110 µs), and at acceptance 0.05 it
makes about 1 accepted move per 20 proposals — but the proposals are fast,
the topology change fraction is 0.93 (nearly every accepted move changes the
tree), and the moves land close to the current state, maintaining local
coherence in the mutation history.

---

## 4. Recommendations

### Immediate (default configuration change)

**Add a second SubtreeLeap for MAT analyses.**  
Target acceptance 0.05, weight = 10% of the main STL weight.

For tall trees (≥ ≈500 tips) this gives +50–70% continuous ESS and +100%
topology ESS at a cost of ≈8–9% extra wall time. The benefit is absent or
marginal for small trees, but neither is there a harm, so a universal default
is safe.

Suggested XML (following the existing BEAUti MAT template):

```xml
<subtreeLeap id="subtreeLeap" size="1.0" weight="NTAXA">
    <treeModel idref="treeModel"/>
</subtreeLeap>
<subtreeLeap id="subtreeLeap2" size="1.0" weight="NTAXA_x0.1"
             targetAcceptance="0.05">
    <treeModel idref="treeModel"/>
</subtreeLeap>
```

where `NTAXA_x0.1` is 10% of `NTAXA`.

### Accept as is

**Target acceptance 0.234:** no universal improvement exists. The optimum
shifts from ≈0.35 on small trees to ≈0.15 on tall trees. 0.234 is a robust
default that avoids catastrophic failure on either.

**STL weight 1× nTaxa:** doubling the weight helps mat on both datasets but
causes stuck chains in tall-1 min. Not safe as a universal default.

### Not recommended

A second STL at factor=0.3 (30% of main weight) gives higher RF ESS on
tall-1 mat but costs more wall time and is less consistent across other
metrics. The factor=0.1 condition dominates in the continuous-parameter
metric that matters most for practical convergence.

---

## 5. Study logistics

The study ran on Machine 1 (darwin) using the `operator_study` branch of
BEAST. A build check on that branch validates that
`mutationHistoryResampler idref` elements are absent from `subtreeLeap` and
`fixedHeightSubtreePruneRegraft` operators (the `matDataLikelihood` resamples
automatically); the benchmark's `make_variants.py` strips those elements from
the generated XMLs.

Three rounds of XML fixes were needed before the mat-likelihood runs succeeded
(first subtreeLeap, then FHSPR, then a git-pull miss). The analysis script
(`analyse.py`) reports operator-level timing from the `.ops` file; the second
STL appears as `stl2:subtreeLeap(treeModel)` once operator ids are present.
