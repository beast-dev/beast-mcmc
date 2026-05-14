# Beanchmark-small
Small benchmark datasets.
< 200 taxa, < 300 sites patterns

## benchmark-small-1.xml

Source data: YFV.fasta
Taxa: 71
Sites: 654
Site patterns: 280

Model: JC69, strict clock, coalescent constant population tree prior
Priors: rate ~ ApproximateReferencePrior, popSize ~ Gamma(0.001, 1000)
Operators: 
  - Scale on rate with weight 3.0
  - Scale on popSize with weight 3.0
  - SubtreeLeap on tree with weight 71.0
  - fixedHeightSPR on tree with weight 7.1
  - upDown on rate and popSize with weight 3.0

MCMC: 10,000,000 states, logging every 10,000 states

## benchmark-small-2.xml

Source data: InfluenzaA_NewYork_HA_2000-2003.nex
Taxa: 165
Sites: 1698
Site patterns: 258

Model: JC69, strict clock, coalescent constant population tree prior
Priors: rate ~ ApproximateReferencePrior, popSize ~ Gamma(0.001, 1000)
Operators:
- Scale on rate with weight 3.0
- Scale on popSize with weight 3.0
- SubtreeLeap on tree with weight 71.0
- fixedHeightSPR on tree with weight 7.1
- upDown on rate and popSize with weight 3.0

MCMC: 10,000,000 states, logging every 10,000 states

