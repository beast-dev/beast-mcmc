## benchmark-big-1.xml

Source data: SARS-CoV-2_B.1.1.7_UK.fasta
Taxa: 976
Sites: 29409
Site patterns: 2078

Model: JC69, strict clock, coalescent constant population tree prior
Priors: rate ~ ApproximateReferencePrior, popSize ~ Gamma(0.001, 1000)
Operators: 
  - Scale on rate with weight 3.0
  - Scale on popSize with weight 3.0
  - SubtreeLeap on tree with weight 976.0
  - fixedHeightSPR on tree with weight 97.6
  - upDown on rate and popSize with weight 3.0

MCMC: 10,000,000 states, logging every 10,000 states

## benchmark-big-2.xml

EBOV_1610_genomes.fasta	
Taxa: 1610	
Sites: 18996	
Site patterns: 7927
Model: JC69, strict clock, coalescent constant population tree prior
Priors: rate ~ ApproximateReferencePrior, popSize ~ Gamma(0.001, 1000)
Operators:
- Scale on rate with weight 3.0
- Scale on popSize with weight 3.0
- SubtreeLeap on tree with weight 1610.0
- fixedHeightSPR on tree with weight 161.0
- upDown on rate and popSize with weight 3.0

MCMC: 1,000,000 states, logging every 1000 states
