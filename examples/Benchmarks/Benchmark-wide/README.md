# Beanchmark-wide
Few taxa, many site patterns benchmark datasets.
< 100 taxa, ~2000 sites patterns

## benchmark-wide-1.xml

Source data: VARV.fasta
Taxa: 64
Sites: 185578
Site patterns: 1955

Model: JC69, strict clock, coalescent constant population tree prior
Priors: rate ~ ApproximateReferencePrior, popSize ~ Gamma(0.001, 1000)
Operators: 
  - Scale on rate with weight 3.0
  - Scale on popSize with weight 3.0
  - SubtreeLeap on tree with weight 64.0
  - fixedHeightSPR on tree with weight 6.4
  - upDown on rate and popSize with weight 3.0

MCMC: 10,000,000 states, logging every 10,000 states


