# Integrated coalescent prior test XMLs (PR #1075)

Paired "full" vs "marginalised" BEAST XMLs reproducing the benchmark scenarios
from beast-dev/beast-mcmc#1075. "Full" places an inverse-Gamma(alpha=2, beta=10)
prior on the constant population size `Ne` (`E[Ne] = 10`) and samples it with a
`scaleOperator`, using the ordinary `coalescentLikelihood`. "Marginalised"
integrates `Ne` out analytically using the new `coalescentLikelihoodIG`
element (`IGCoalescentLikelihood`/`IGCoalescentLikelihoodParser`), with no
population-size parameter or operator at all.

## Scenarios

| Files                                    | Taxa | Sequence data | Notes |
|-------------------------------------------|------|----------------|-------|
| `prior_n5_{full,marginalised}.xml`        | 5, contemporaneous  | none | prior-only; theoretical `E[root height] = 2*beta*(n-1)/((alpha-1)*n) = 16` |
| `prior_n50_{full,marginalised}.xml`       | 50, contemporaneous | none | prior-only |
| `prior_n17_hetero_{full,marginalised}.xml`| 17, heterochronous  | none | prior-only; same tip dates as `examples/TestXML/TreePriors/testConstantSize.xml` |
| `posterior_n17_hetero_{full,marginalised}.xml` | 17, heterochronous | real Dengue-4 alignment | full posterior; adapted from `examples/TestXML/TreePriors/testConstantSize.xml` |

The `prior_*` XMLs regenerate from `generate_prior_xmls.py`; the
`posterior_n17_hetero_*` pair was hand-adapted from the real 17-taxon Dengue-4
example already in this repository (`testConstantSize.xml`) by swapping only
the population-size/coalescent-prior block — everything else (alignment,
substitution model, clock, operators) is unchanged between the two.

## Running

```sh
java -cp <beast-classpath> dr.app.beast.BeastMain -overwrite prior_n5_full.xml
```

The `posterior_n17_hetero_*` pair uses `treeDataLikelihood`, which normally
needs the native BEAGLE library. Without it installed, force the pure-Java
fallback:

```sh
java -Djava.only=true -cp <beast-classpath> dr.app.beast.BeastMain -overwrite posterior_n17_hetero_full.xml
```

All 8 files were verified to parse and run to completion (a few iterations)
against this branch's build before committing.

Chain length is set to 10,000,000 (matching the PR), logging every 1,000
states. Adjust `chainLength`/`logEvery` in the `<mcmc>` element for shorter
smoke tests.
