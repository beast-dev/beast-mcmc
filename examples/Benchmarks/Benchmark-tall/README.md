# Beanchmark-tall
Many taxa, few site patterns benchmark datasets.

## benchmark-tall-1.xml

Source data: H3N2_HA.fasta
Taxa: 1441
Sites: 987
Site patterns: 593

Model: JC69, strict clock, coalescent constant population tree prior
Priors: rate ~ ApproximateReferencePrior, popSize ~ Gamma(0.001, 1000)
Operators:
- Scale on rate with weight 3.0
- Scale on popSize with weight 3.0
- SubtreeLeap on tree with weight 1441.0
- fixedHeightSPR on tree with weight 144.1
- upDown on rate and popSize with weight 3.0

MCMC: 10,000,000 states, logging every 10,000 states

