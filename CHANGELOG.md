# BEAST X Version History

Bayesian Evolutionary Analysis Sampling Trees  
https://beast.community/

All issues are tracked at https://github.com/beast-dev/beast-mcmc/issues

---

## BEAST X v10.6.0

Changes since the v10.5.0 release (2 July 2025).

### Compatibility notes

- **`<birthDeathModel>` has been renamed `<Gernhard08BirthDeathModel>`** as part of
  the split of the birth-death models into their own package (`dr.evomodel.birthdeath`).
  The old element name is not registered as an alias, so XML written by earlier BEAUti
  versions with a birth-death tree prior will fail to parse until the element is renamed.
  `<newBDSSHistorySimulator>` is likewise now `<ebdsHistorySimulator>` with no alias.
  `<newBirthDeathSerialSampling>` is renamed `<episodicBirthDeathSamplingModel>` but the
  old name is kept as an alias.
- **Checkpoint files** written by earlier versions cannot be expected to resume a run
  exactly: tips are now numbered by taxon name rather than by position in the starting
  tree ([#1225](https://github.com/beast-dev/beast-mcmc/issues/1225)), operator state is
  matched by operator name rather than by position, and adaptive multivariate normal
  operators now require their saved `operator_state` lines. A run started with the same
  seed as under v10.5.0 follows a different trajectory.
- **`<structuredCoalescent>`** (BASTA) now requires the population sizes inside a
  `<popSizes>` wrapper element rather than as a bare `<parameter>` child, and
  `useAmbiguities` defaults to `true`.
- **`<timeVaryingRates>`** now requires an `<IntegratedTransformedSplines>` child; the
  former `coefficients`/`knots`/`degree` content is no longer accepted.
- **`<freeRateSimplexTransform>`** is the new name of `<realDifferencesToSimplexTransform>`.
- **Taxon names** that need quoting in tree and NEXUS files are written with double
  quotes rather than single quotes, and a taxon id containing a double quote is rejected.
- **Numerical gradients** (used by HMC gradient checks and the numeric fall-backs) use a
  central-difference step of ε^(1/3)·max(|x|, 1) rather than √ε·(|x|+1), so reported
  numeric gradients differ slightly from earlier versions.
- **`Q.3Di.FSlogodds`** is the new name of the `ThreeDi` amino acid model (the old name
  is accepted); its matrix was regenerated from FoldSeek's log-odds matrix and the values
  differ from v10.5.0.
- **The FreeRate site model** logs its rate and weight columns sorted by increasing rate
  rather than by weight, and prefixes the column names with the model's id
  (`siteModel.rates1`) so that several site models can be logged together.

### New features

#### BEAST: new models

- **Q.3Di.AF and Q.3Di.LLM substitution models** for FoldSeek's 3Di structural alphabet,
  selected with `type="Q.3Di.AF"` or `type="Q.3Di.LLM"` in `<empiricalAminoAcidModel>`.
  Each carries its own equilibrium frequencies (Garg and Hochberg, 2025, Mol Biol Evol
  42:msaf124). A `<frequencies>` child of `<empiricalAminoAcidModel>` is now honoured;
  it was silently ignored before.
- **FreeRate among-site rate variation.** `<freeRateSiteRateModel>` takes an optional
  `<proportionInvariant>` parameter (a zero-rate category is added in front of the free
  categories) and its `<substitutionModel>` child is now optional. The model is loggable,
  giving the sorted rates, weights and a `rates_x_weights` column. New operators
  `<freeRateUpDownOperator scaleFactor=... weight=...>` (scales two rates keeping the
  weighted mean at 1) and `<freeRateDeltaExchange delta=... weight=...>` (a logit-scale
  move on a pair of weights with the compensating rate change), both taking `<weights>`
  and `<rates>` children; `<freeRateSimplexTransform dimension=...>` maps positive reals
  to the weighted simplex; `<newTransformedParameterOperator>` wraps any operator with the
  Jacobian correction for a transform. Transformed multivariate parameters are
  checkpointed by their untransformed dimension. The two FreeRate operators are
  registered in `release_parsers.properties`, and `<freeRateSiteRateModel>` requires its
  `categories` attribute and declares the `<substitutionModel>` and
  `<proportionInvariant>` children, so XML using them passes `-strict` parsing.
- **Integrated coalescent prior.** `<coalescentLikelihoodIG>` marginalises the constant
  population size under an inverse-gamma(`alpha`, `beta`) prior, with `<include>`/`<exclude>`
  taxon lists; ported from Luiz Max Carvalho's work
  ([#1075](https://github.com/beast-dev/beast-mcmc/pull/1075)) and corrected for a
  population size shared across trees
  ([#1273](https://github.com/beast-dev/beast-mcmc/pull/1273)).
- **Gaussian process coalescent.** `<multiLocusNPCoalescentLikelihood>` is the coalescent
  density alone, separated from the GMRF prior that used to be built into
  `GMRFMultilocusSkyrideLikelihood`, so the log population sizes can be given any prior.
  It takes one or more `<populationTree>` elements, `<populationSizes>` (optionally with
  a transform applied to the starting values), either `<gridPoints>` or `<cutOff>` plus
  `<numGridPoints>`, and an optional per-tree `<ploidy>`.
  `<multilocusNPCoalescentLikelihoodGradient>` gives its gradient for HMC and
  `<popSizesLogger>` logs population sizes on the natural scale (`PopSize.1..N`) for this
  or the GMRF skyride and skygrid likelihoods. The
  `<gaussianProcessField>` prior now takes standalone `<basis>` elements (a design matrix
  or weight provider plus a `<kernel>`, an optional `<weightFunction>` of type `sigmoid`,
  `sigmoidComplement`, `identity` or `linear` for non-stationary kernels, and
  `orthogonalProjection` / `unitaryVariance` attributes), a `NeuralNet` kernel, an
  estimable nugget, `<gaussianProcessKernelGradient>` for the `scale` and `length`
  hyperparameters, and `<gaussianProcessConditionalDerivative>` for the derivative of the
  field conditional on its values. `<gaussianProcessBasisApproximation>` is a
  reduced-rank Hilbert-space approximation (`degree` selects the exponential, Matérn 3/2,
  Matérn 5/2 or squared-exponential spectral density) and `<BSplines>` a B-spline prior.
  `GMRFMultilocusSkyrideLikelihood` reports its coalescent, field and total
  log-likelihoods.
- **New clock models.** `<nonParametricRates>` integrates a squared Gaussian process
  basis approximation over each branch; `<nonparametricRates2>` does the same with
  `<IntegratedTransformedSplines>`; `<gridBasedBranchRateModel>` takes `<gridPoints>` and
  `<levelSpecificRates>` (one more rate than grid points) with
  `<gridBasedBranchRateModelGradient>`; `<logLinearBranchRate>` averages
  exp(`slope`·t + `intercept`) over each branch; `<timeVaryingRates>` gains the
  `piecewiseLogLinear` functional form (with an optional `<slope>` for the last epoch) and
  the `Integrable` form driven by `<IntegratedTransformedSplines>` (children
  `<coefficients>` and `<intercept>`, attributes `knots`, `degree="2|3"` and
  `transform="squared|exponential"`).
- **Strongly lumpable CTMC models.** `<stronglyLumpableCtmcRates>` builds the rates of a
  discrete-state CTMC whose states are partitioned into `<lump>` elements (each a
  `<stateSet>`, within-lump `<rates>` and `<proportions>` towards other lumps), with the
  across-lump `<rates>` at the top level and an optional `normalize` attribute. It is
  used through `<logRateSubstitutionModel><rateProvider>`; `<approximateLogCtmcRateGradient>`
  returns the lumpable gradient for HMC, and `<stronglyLumpableGenerator>` writes a
  template XML block (`generated_strongly_lumpable_block.xml`) from a data type and
  state sets.
- **Transformed-rate substitution models.** `<logRateSubstitutionModel>` (alias
  `<transformedRateSubstitutionModel>`) accepts `<transformedRates>` (a parameter plus a
  transform) or a `<rateProvider>` instead of `<logRates>`, implements exact
  differentiation with respect to each rate, and `<approximateLogCtmcRateGradient>` has
  the alias `<approximateTransformedCtmcRateGradient>`. `<complexSubstitutionModel>`,
  `<glmSubstitutionModel>` and `<logRateSubstitutionModel>` take
  `scaleRatesByFrequencies="false"` to build Q without multiplying by the equilibrium
  frequencies. A BSSVS indicator with a value other than 0 or 1 is now a parse error.
  `<infinitesimalRatesLogger>` takes `diagonalElements`, `order` and a transform.
- **Time-varying uncertain tip states.** `<timeVaryingFrequencies>` places a prior on
  the states of an unobserved tip using the frequency vector of the epoch (`<cutOff
  time=...>`) its height falls into, and `<tipStateOperator>` resamples those states
  (Suchard and Pekar, in preparation).
- **Markov-modulated models** can be used in branch-specific substitution-parameter
  branch models and provide a first-order approximate gradient;
  `<branchSubstitutionParameterGradient>` passes `mode` through to the branch-specific
  case and `dim="-1"` differentiates all dimensions; exact derivatives are correct for
  complex eigenvalues.
- **Continuous traits.** `<taxonEffectTraitDataModel>` adds an estimable additive effect
  per taxon (`<effects>`, optionally shared with another model, with a single tree-wide
  `<sign>`), with `<taxonEffectGradient>` and support in `<varianceProportionStatistic>`;
  `<sliceData by="rows|columns">` exposes the rows or columns of a matrix parameter as
  parameters; `<bandedMultivariateDiffusionModel representation="dense|sparse">` builds a
  banded precision from a `<precisionMatrix>` block replicated `replicates` times or
  structured by a GMRF, for very high-dimensional traits, with sparse branch-rate and tip
  gradients and `<blockDiagonalMatrixParameter>`; `<sampledCategoricalVarianceModel>`
  gives each trait a category-specific variance with its gradient;
  `<traitDataLikelihood forceCompletelyObserved="true">` treats the data as fully
  observed; the correlation matrix of a factor model is available as a parameter
  (`<covarianceToCorrelation>`).
- **Structured coalescent (BASTA).** The `gradientImplementation="adjoint|legacy"`
  attribute of `<structuredCoalescent>` selects the new adjoint (reverse-mode) delegate;
  `<structuredCoalescentLikelihoodGradient wrtParameter=...>` supports `migrationRate`
  for any `GeneralSubstitutionModel` (with BSSVS indicators applied through the chain
  rule), `populationSize`, `glmCoefficients` and `glmRandomEffects`. Per-deme
  `<growthRates>` give exponential population growth; population size dimensions are
  named by state (`popSizes.Belgium`). Ancestral state reconstruction, sampling of
  ambiguous tip states, and `<counts>`/`<rewards>` Markov jump and reward registers
  (with `<markovJumpsLogger name=... reportUnconditional=...>`) are supported, with the
  old reconstruction kept as the `<stateTagName>_old` trait. Trees with travel-history
  (ancestral trait) nodes are accepted. The native library is now `hmsbeagle-bit`;
  `-beagle_basta_threads N` / `-beagle_basta_threading_off` control its CPU threading
  and `-beagle_threading none|cpp|openmp` the general BEAGLE threading model.
  Citation: Shao et al. (2026) PNAS 123:e2602412123.
- **Birth-death models.** `<birthDeathLikelihood>` and `<birthDeathLikelihoodGradient>`
  (release parsers) take a `<model>` and `<tree>`, with `<include>`/`<exclude>` and
  `useNewLoop` for the efficient implementation. `<episodicBirthDeathSamplingModel>`
  can be parameterised by `<R0>`, `<D>` and `<S>` instead of birth, death and sampling
  rates; `treatmentProbability` and `samplingProbability` are now required.
- **Finite mixture parameters.** `<finiteMixtureParameter>` maps `<values>` onto
  `<categories>` and `<finiteMixtureModel>` (attributes `prior="geometric|poisson"`,
  `hyperParameter` and `constraint="none|ordered"`) places a prior on the number of
  occupied clusters, logging
  sorted values and categories (used to cluster skyline population sizes).
- **Recombination height bound statistic.** `<recombinationHeightBoundStatistic>` takes
  two or more `<recombinant>` taxon sets and reports the maximum parent height and the
  minimum MRCA height across them.
- **Distributions and transforms.** `<cauchyDistributionModel median scale>` with a
  gradient (and `<cauchyPrior>` rewritten with closed forms and a gradient);
  `<affineTransform location scale dim>` for (x − location)/scale;
  `<UnitSimplexTransform dimension>` (stick-breaking) and `<UnitSimplexViaSoftmaxTransform>`;
  `type="sigmoid|squared|abs|affine"` in `<transform>`; composition of two multivariable
  transforms.
- **Laplace marginal likelihood.** `<approximateTreeDataLikelihood>` wraps a
  `<maximizeWrtParameter>` and returns the Laplace approximation to the log marginal
  likelihood from the diagonal Hessian at the maximum.
- **Geodesic HMC.** `<geodesicHamiltonianMonteCarloOperator>` samples on spheres and
  Stiefel manifolds, with masks.
- **Coalescent utilities.** `<coalescentRescaler>` rescales a tree simulated under a
  constant population to a piecewise-constant demographic; `<piecewisePopulation>` accepts
  `<gridPoints>` in place of epoch widths, an optional transform on `<populationSize>`,
  and omits `<growthRates>` for the piecewise-constant case; `<ctmcScalePrior>` takes
  several trees and a `proportion` of their mean length.

#### BEAST: run-time features

- **Output directory** — `-output <directory>` on the command line, or "Write all
  output files to a folder" in the launch dialog, writes every log, tree, checkpoint and
  citation file to the chosen directory, checked to exist and be writable before the run
  starts. A checkpoint file named in the XML is written next to the log files.
  ([#1267](https://github.com/beast-dev/beast-mcmc/pull/1267)) ([#1280](https://github.com/beast-dev/beast-mcmc/pull/1280))
- **Launch dialog** — the BEAGLE options are shown only when the BEAGLE library is found,
  otherwise the dialog explains where to get it. Quit and the window's close button quit
  the application, and Run without an XML file returns to the dialog.
- **Fast tree model** — `<treeModel implementation="fast">` keeps the tree in flat arrays
  stored copy-on-write, producing the same MCMC output as the default tree model at a
  lower cost per step on large trees. Node rates and node traits are not supported.
  ([#1274](https://github.com/beast-dev/beast-mcmc/pull/1274))
- **Faster coalescent calculations** — coalescent interval times are sorted incrementally
  rather than from scratch at each evaluation, and the CTMC scale prior keeps its tree
  length up to date incrementally. ([#1275](https://github.com/beast-dev/beast-mcmc/pull/1275))
- **Operator analysis and logging** — the operator analysis is written as CSV when the
  operator log file name ends in `.csv` or `format="csv"` is set, and includes the
  weights; the operators element can be logged (`_time`, `_size` columns); tree operators
  record the size of each move and whether the topology changed (`topologyChanged`,
  `mean_delta` columns). Operator profiling is in nanoseconds.
- **Chain profiling** — `-Dmcmc.profile=steps`, `likelihoods` or `all` times the phases
  of every MCMC step per operator and every likelihood evaluation, reported after the
  operator analysis. ([#1276](https://github.com/beast-dev/beast-mcmc/pull/1276))
- **The command line help is grouped into sections** — the options are listed under
  headings (General, The chain, Site patterns, BEAGLE, Running several chains,
  Checkpointing, Citations and plugins, Information) rather than as one long list.
- **Overriding the chain length and logging frequency** — `-chain_length` (short `-cl`)
  runs for a given number of states and `-log_every` (short `-le`) logs every given number
  of states, whatever the file says, so that a file can be run short, or sampled more or
  less often, without being edited. Both are reported as the run starts. The logging
  frequency applies to every logger that has one, the screen log included, and a logger
  without a fixed frequency is left alone. The options set the `mcmc.chain_length` and
  `mcmc.log_every` properties, which can also be given directly with `-D`.
- **Checkpointing** — adaptive multivariate normal operators save and restore their full
  adaptive state; runs with an `<empiricalTreeDistributionModel>` save the current tree
  index; `-checkpoint_seed <seed>` sets the seed used when a checkpoint is loaded (the
  saved random number generator state is restored otherwise); `<prngStateLogger>` logs
  the generator state for debugging.
- **BEAGLE** — `-beagle_tensor_core` requests tensor cores on a GPU;
  `-beagle_threading none|cpp|openmp` selects the threading implementation.
- **Numerical gradients** — `<parallelNumericGradient>` evaluates a numeric gradient over
  replicate likelihoods on a thread pool; `<numericalGradient>` has the alias
  `<numericGradient>`; gradient checks (including `<skylineGradient>` and
  `<skyGlideGradient>`) take `numericStepSize="tiny|small|medium|large"` or a number.
- **Empirical tree distributions** work with the fast coalescent intervals and the
  efficient speciation likelihood.
- **Epidemiological statistics** — `<epidemiologyStatistic>` can be written as
  `<doublingTime>`, `<growthRate>` or `<R0>`, taking a `<growthRate timeUnits=...>` or
  `<doublingTime timeUnits=...>` parameter and, for R0, a `<serialInterval mean stdev>`.
- `duplicatedParameter` can take a fixed, constant number of copies as an attribute;
  `glmModel` takes an id used when logging its rates, and GLM coefficient and indicator
  dimensions are named after their predictors.
- Every command-line tool built on the shared `Arguments` class accepts short
  single-dash option forms and lists both forms in its help.
- On macOS, BEAST, TreeAnnotator and LogCombiner have the application menu (About, Quit)
  and the Dock icon, and BEAST's About box is laid out like BEAUti's.

#### BEAUti

- **Free Rates site heterogeneity model** — the Sites panel's "Site Heterogeneity Model"
  offers "Free Rates" for nucleotide and amino acid partitions, with 2 to 16 rate
  categories (the "Number of Rate Categories" menu offers 4 to 16 for the gamma model,
  as before). The XML has a `<freeRateSiteRateModel>` with `freeRates.rates` and
  `freeRates.weights` parameters, a `freeRates.weightedRates` product statistic that
  carries the Dirichlet(1) prior on the rates, a symmetric Dirichlet(2) prior on the
  weights (the concentration above one keeps a category's weight from collapsing to
  zero, where its rate is unidentifiable and the likelihood numerically unstable), the
  `<freeRateDeltaExchange>` and `<freeRateUpDownOperator>` operators, and the site
  model in the log so that the rates and weights are logged sorted by rate. The model
  can be linked or unlinked across codon positions like the gamma model. Generalized
  stepping-stone sampling has no working prior for it yet and is refused.
- **Dirichlet prior concentration** — the prior dialog has a "Concentration" field for
  every Dirichlet prior (base frequencies, GTR rates, the free rates weights and
  rate shares), written as the `alpha` of `<dirichletPrior>`; it was fixed at 1.0.
- **Structured coalescent (BASTA)** — a discrete trait partition's "Discrete Trait
  Substitution Model" can be "Backward-in-time CTMC", with a constant-population
  coalescent model, symmetric, asymmetric or GLM structure, BSSVS, and a "Share coalescent
  model across demes" option. The tree prior is then set by the structured coalescent.
  The generated XML writes `<structuredCoalescent type="BASTA">` with a
  `<complexSubstitutionModel normalized="true" scaleRatesByFrequencies="false">`,
  `backward.<model>.rates` and `.indicators` parameters, a gamma prior on
  `structuredCoalescent.popSizes` (shared through `<duplicatedParameter>` or per deme),
  and the matching operators and log columns. Root frequencies and state-change counting
  are not available for the backward-in-time model.
- **Fixed local clock** — the Clocks panel shows the Yoder & Yang (2000) citation and
  lists the taxon sets that will be given their own rate (those marked monophyletic), or
  explains how to define one. The "Subtree prior model" drop-down is greyed out when no
  monophyletic taxon set exists.
- **Epidemiological statistics** — for the exponential, logistic and expansion coalescent
  priors, options to log the doubling time from the growth rate (or the growth rate from
  the doubling time) and R0 from the growth rate given a gamma-distributed serial
  interval mean and standard deviation.
- **Gao et al. (2023) prior specification for discrete trait models** — an option in the
  Sites panel replaces the default priors on a forward-in-time discrete trait model with
  the more diffuse specification of Gao, May, Rannala & Moore (2023, PNAS
  120:e2213913120), including a new "Exponential with estimated mean" prior on the
  overall rate scalar.
- **macOS integration** — the application menu, the Dock icon and files opened from the
  Finder go through `java.awt.Desktop`, the title bar takes the window's colour, the
  toolbar icons are drawn as vectors for high-resolution displays, combo boxes use Aqua's
  pop-up style and highlight the item under the mouse, and `-Dswing.defaultlaf` overrides
  the system look and feel.
  ([#1277](https://github.com/beast-dev/beast-mcmc/pull/1277)) ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))
- Import errors are reported in a dialog that says what is wrong with the file, and the
  macOS file dialog accepts several files at once.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))
- The operator analysis log is off by default.
- Birth-death tree priors are written as `<Gernhard08BirthDeathModel>`.

#### TreeAnnotator

- **Multiple HPD intervals**: `-hpdIntervals 0.95,0.90,0.5` (short `-hi`), with
  `-hpdLimit` the minimum number of height values for an interval (default 10).
- **Kernel density estimates of node heights** for visualisation in PearTree: `-kde`,
  with `-kdeCount` curve points (default 100) and `-kdeLimit` minimum height values
  (default 200). Off by default in the GUI until it has an interface.
- **Common ancestor heights**: `-mrcaHeights` (short `-ca`) computes clade MRCA heights
  with the method of Heled and Bouckaert (2013) to avoid negative branch lengths.
- Reports the number of clades above the credibility threshold, the median clade
  credibility of the sampled trees, and the log clade credibility of a user target tree.
- Command-line changes: every option has a short form; `-limit` is now
  `-limitFrequency` (default 0.0) and `-limitCount` defaults to 0 rather than 5;
  `-heights keep` and the `mrc` tree type are removed; the annotated tree is written
  without the NEXUS taxa block unless `-nexusHeader` is given; `-forceDiscrete` and
  `-ess` are withdrawn. The HIPSTR citation is printed.

#### Other tools

- **LogCombiner burn-in** — the burn-in can be given as a number of trees or log samples
  (`-burninTrees`, short `-bt`) or as a percentage of each file (`-burninPercent`, `-bp`)
  as well as in MCMC states (`-burnin`); the dialog has a "Burn-in given as" choice and the
  table column names the unit. For every file LogCombiner now reports how many trees or
  samples were read and written, and warns when the burn-in removed nothing, only the first
  sample or everything, which is what happens when a number of trees is entered where states
  are expected. Lines that could not be read, which are skipped as a possible corruption,
  are counted and reported separately rather than being blamed on the burn-in.
  ([#1247](https://github.com/beast-dev/beast-mcmc/issues/1247))
- TreeStateTimeSummarizer writes `startHeight`, `endHeight` and `branchingEvents`
  columns for each state interval.
- `scripts/fit_rate_prior.py` fits a lognormal distribution to a set of HPD bounds for
  use as a rate prior.
- SequenceSimulator handles branch-specific substitution models and the current site
  rate model classes.

### Bug fixes

#### BEAST

- **Fixed a run resumed from a checkpoint not following the run that wrote it** — tips
  were numbered by their position in the (random) starting tree, so a resumed run had
  different node numbers and a permuted `branchRates.categories`, and adopting the saved
  tree permuted the categories again. Tips are now numbered by taxon name, in both the
  default and the fast tree model, and the saved tree is adopted without permuting
  node-indexed parameters. ([#1225](https://github.com/beast-dev/beast-mcmc/issues/1225))
- **Fixed the branch-rate diagonal Hessian for transformed rates** — the reciprocal, exp
  and location-scale transforms gave a wrong mass matrix to HMC with (adaptive) diagonal
  preconditioning. Reporting a `traitDataLikelihood` with `arbitraryBranchRates` no longer
  fails with an ArrayIndexOutOfBoundsException.
  ([#1279](https://github.com/beast-dev/beast-mcmc/pull/1279))
- **Fixed site pattern compression taking quadratic time** — a 70 taxon by 644,000 site
  alignment never finished loading in BEAST or BEAUti. Patterns are found by a hash lookup.
  ([#1269](https://github.com/beast-dev/beast-mcmc/pull/1269))
- **Fixed BASTA giving -Inf for every state on an ancestral-trait tree**, where a
  travel-history ancestor sits on a zero-length coalescent interval by construction.
- **Fixed BigFastTreeIntervals** restoring stale intervals after a rejected proposal,
  giving wrong intervals for tip-height moves after a full rebuild, and ordering the root
  before other nodes of the same height.
- **Fixed the coalescent gradient** miscounting zero-length intervals when internal and
  external nodes share a height, and the skyglide gradient with respect to the log
  population size when the slope is near zero.
- **Fixed the smooth skygrid likelihood** when the root lies beyond the last grid point.
- **Fixed the GMRF** scale matrix omitting the precision when the field determinant is
  factored, not listening to its weight provider, and the determinant for proper fields,
  which now uses a tridiagonal recursion.
- **Fixed pre-order partials indexing in the BEAGLE gradient delegates**, which gave
  wrong results after a restore with double buffering.
- **Fixed Markov jumps not being redrawn after a BSSVS indicator change.**
- **Fixed `arbitraryBranchRates` randomisation** under the reciprocal transform
  returning −exp(x) rather than exp(−x).
- **Fixed FreeRate with a proportion of invariant sites**, the FreeRate delta exchange
  Hastings ratio, and FreeRate checkpointing.
- **Fixed `<ancestralTreeLikelihood>`, `<markovJumpsTreeLikelihood>` and
  `<treeLikelihood>` requiring a `GammaSiteRateModel`**, so a `<siteModel>` without gamma
  categories (which produces a `DiscretizedSiteRateModel`) could not be used, and
  `<treeTrait>` now reports a child that is not a trait provider.
- **Fixed the `dNdSLogger` for prefixed robust counting partitions.**
- **Fixed `<empiricalDistributionLikelihood>` parsing for radiocarbon dating.**
- **Fixed `<R0>` in `<epidemiologyStatistic>` computing a doubling time**, the gamma
  shape of the serial interval, and the logged column name.
- **Fixed InverseGammaDistribution counting its normalising constant twice** in `pdf`
  and `logPdf`.
- **Fixed the latent liability model's variable sign** being taken from the starting
  value, which threw when a variable started exactly on zero; it is taken from the bounds.
- **Fixed a NullPointerException in the rate covariance statistic** with the fast tree
  model.
- **Fixed LogCombiner taking the wrong samples when resampling at a lower frequency** —
  the samples to keep were chosen by a running count of the samples kept so far rather than
  by their own state number, so a burn-in shifted the choice by one sampling interval and
  the samples taken from each file depended on how many had been taken from the files before
  it. Combining with a burn-in of 1,000,000 and a resampling frequency of 100,000 took the
  sample at state 1,090,000 rather than 1,100,000. The state numbers written were that same
  running count, so with a burn-in they did not correspond to the samples beneath them.
  Samples are now chosen by their own state number relative to the first one kept from the
  same file, and the state numbers written run consecutively through the combined output.
  Combining without resampling is unaffected.
  ([#857](https://github.com/beast-dev/beast-mcmc/issues/857))
- **Fixed LogCombiner failing with an ArithmeticException** when renumbering output states
  without also resampling.
- **Fixed LogCombiner dropping the first sample of the first log file** (the state 0 line
  was compared against a column count that was only set from the second line), and having
  two options with the short name `-r`.

#### BEAUti

- **Fixed BEAUti reverting clock-rate prior to CTMC reference prior** — custom rate
  priors (e.g. a lognormal prior on the substitution rate) were silently replaced by the
  CTMC reference prior when the XML was written.
  ([#1238](https://github.com/beast-dev/beast-mcmc/issues/1238))
- **Fixed lognormal prior parameters not being saved** from the prior dialog.
  ([#1258](https://github.com/beast-dev/beast-mcmc/issues/1258))
- **Fixed log and tree file names being truncated** when the name stem was too long.
  ([#1259](https://github.com/beast-dev/beast-mcmc/issues/1259))
- **Fixed saving and opening a BEAUti document** — saving failed for any analysis with an
  imported tree, and a document that did open had its data types replaced by copies so
  the Partitions and Site Models panels still showed the previous document.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))
- **Fixed hierarchical models surviving a change to the analysis** — a hierarchical
  phylogenetic model holds its parameters directly, and changing the analysis, such as
  partitioning into codon positions differently, does not destroy the parameters it
  replaces, it only stops them being used. A model built beforehand was still written out,
  so the file referred to parameters that it did not declare and BEAST stopped with
  "Object with idref=CP1.kappa has not been previously declared". Such a model now has its
  departed parameters dropped before the file is written, and is removed altogether if
  fewer than two are left; a parameter still in the analysis that is left without a model
  is given an ordinary prior, and the user is told what was dropped and taken to the Priors
  panel. ([#1062](https://github.com/beast-dev/beast-mcmc/issues/1062))

- **Fixed a NEXUS file with character sets that do not fit its alignment failing to
  import** — a character set running past the end of the alignment, starting before its
  first site, or ending before it starts, failed while the sites were being taken, with an
  index error that said nothing about which set was at fault; in the application the import
  simply did not happen. Such a file is now refused with a message naming the character
  set, the range it asks for and the length of the alignment.
  ([#1066](https://github.com/beast-dev/beast-mcmc/issues/1066))

- **Fixed the fixed local clock failing to generate XML** for a taxon set imported from a
  one-column file. ([#1231](https://github.com/beast-dev/beast-mcmc/issues/1231))
- **Fixed every taxon set being nested as a subtree of the random starting tree**, which
  made BEAST fail with "subtrees' taxa overlap" when sets overlapped. Only taxon sets
  marked monophyletic are nested.
- **Fixed per-tree taxon lists not being written when partitions have different taxa.**
  ([#1271](https://github.com/beast-dev/beast-mcmc/pull/1271))
- **Fixed the HMC relaxed clock and shrinkage local clock XML with several clock models.**
  ([#1272](https://github.com/beast-dev/beast-mcmc/pull/1272))
- **Fixed the gamma and Cauchy relaxed random walks** placing the gamma on a scalar of
  the diffusion variance rather than the precision, so the Cauchy setting was not a
  Cauchy. ([#1279](https://github.com/beast-dev/beast-mcmc/pull/1279))
- **Fixed a bit-flip operator on the rate indicators being written without BSSVS.**
  ([#1278](https://github.com/beast-dev/beast-mcmc/pull/1278))
- **Fixed switching a discrete trait model from backward-in-time back to forward-in-time**,
  and the forward/backward XML generation.
- **Fixed `includeStem` always written as true** for taxon set TMRCA statistics.
- **Fixed robust counting XML** for prefixed traits and for alignments with sites that
  are gaps in every sequence.

- **Renaissance counting is no longer offered for a partition that is not a whole number
  of codons long** — the three codon positions are compared site by site, so they must hold
  the same number of sites, which taking every third site only gives when the partition's
  length is a multiple of three. BEAUti offered the option regardless, and BEAST then
  refused the file with "Codon-partitioned robust counting requires all partitions to have
  the same length". The option is now disabled with the reason shown, and a file whose
  option was set before this was checked is reported rather than written.
- **Fixed amino acid site model XML generation** failing with a NullPointerException,
  and marginal likelihood XML for the 3Di models.
- **Fixed importing a traits file turning on tip dates** for the whole analysis.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))
- **Fixed names not being escaped in the XML** — a partition, model, taxon set, trait or
  file name containing an ampersand, an angle bracket or a quote produced a file BEAST
  could not read. ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))
- Generate no longer leaves a truncated XML file behind when it fails, takes the value
  being edited on the current panel before validating, and reports a bad combination of
  models before anything is written, switching to the panel concerned.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))
- Setting the Clock or Tree column of the Partitions table follows the same rules as the
  link buttons, and unlinking or cloning a model copies the complete settings.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))
- Dialogs are disposed when closed, imported files are closed, clearing a number field or
  selecting a separator in the tree prior list no longer throws, hovering over the MCMC
  panel's checkboxes no longer marks the document as changed, and declining to take a
  prior out of a hierarchical model abandons the edit.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))
- The stochastic Dollo model, which produced an XML whose likelihood took no part in the
  MCMC, has been removed. ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))
- Taxon names containing a double quote are rejected on import, with the offending names
  listed, since such names cannot be written to a tree file.

#### TreeAnnotator

- **Fixed TreeAnnotator not summarising 2-dimensional locations correctly** — continuous
  phylogeographic location traits (e.g. `location1` / `location2`) were renamed to
  `location.rate1` / `location.rate2` in the annotated tree.
  ([#1242](https://github.com/beast-dev/beast-mcmc/issues/1242))
- **Fixed performance regression in the "Collecting node information" stage** —
  processing of large tree sets became increasingly slow and eventually stalled.
  ([#1233](https://github.com/beast-dev/beast-mcmc/issues/1233))
- **Fixed continuous phylogeography with fixed trees** — coordinate annotations were
  replaced by branch lengths and discarded when the input topology was fixed.
  ([#1255](https://github.com/beast-dev/beast-mcmc/issues/1255))
- **Fixed a NullPointerException when a user target tree contains clades absent from the
  sampled trees** — such clades get a posterior of 0 and keep their target-tree height.
- Fixed wrong attribute names in the annotated tree, and `-heights keep` being accepted
  for HIPSTR trees, where it is not meaningful.
- **Fixed BEAGLE not being found on Apple Silicon Macs** — BEAST looks in
  `/opt/homebrew/lib` as well as `/usr/local/lib`. Homebrew installs there on arm64, so a
  brew-installed BEAGLE was invisible to both `bin/beast` and BEAST.app, and the failure is
  quiet: BEAST reports no BEAGLE resources and runs the slower Java likelihood rather than
  stopping.

### Other changes

- Updated the HIPSTR and BASTA citations to the published papers
  ([#1243](https://github.com/beast-dev/beast-mcmc/issues/1243)).
- Taxon names that need quoting in tree and NEXUS files (whitespace or any of
  `. ( ) [ ] : ; , = '`) are delimited with double quotes in all of BEAST's output, so
  names containing an apostrophe can be written and read back.
- BEAGLE is loaded with `System.loadLibrary`, which JDK 24 and later warn about unless
  native access is enabled: the `beast.jar` manifest enables it for `java -jar`, and the
  `beast` script adds the option when the installed Java accepts it.
- The warning that a likelihood was created but not used by the MCMC now says "(no ID)"
  for a likelihood that has no id, rather than "(id: null)".
  ([#1019](https://github.com/beast-dev/beast-mcmc/issues/1019))
- Stopping a run uses `Thread.interrupt` rather than the removed `Thread.stop`.
- Development builds report their version as 10.x.x rather than an unexpanded `${version}`.
- **Packaging** — releases are built by a GitHub Actions workflow (`build_release.xml`)
  triggered by a version tag: notarised macOS DMGs (with and without a bundled JVM, built
  with `jpackage`), a Windows `jpackage` image sharing one Temurin runtime between BEAST X,
  BEAUti, TreeAnnotator and LogCombiner, and a Linux tarball that includes the BEASTGen
  shell scripts. The applications are named "BEAST X".
- **Heap sizes** are set as a fraction of the machine's memory (`-XX:MaxRAMPercentage=50`)
  rather than the fixed `-Xmx` values each launcher carried, which ranged from 1024m to
  8192m and differed between the script and the packaged application for the same tool. A
  few thousand sequences of a few tens of thousands of sites exhausted BEAST's old 2 GB
  while the likelihood was still being set up, before the chain started. The percentage is
  a ceiling rather than a reservation, so a small analysis still uses very little, and it
  keeps pace with whatever machine BEAST is run on. `BEAST_XMX=8g beast ...` fixes a size
  instead. On a machine with less than 8 GB, TreeAnnotator now has a lower ceiling than its
  former 4096m.
- The Mac package ships the same `bin/` scripts as the Linux one. The Mac copies were
  identical but for two stale heap sizes, so there is one set to maintain rather than two.
- BEAUti has a headless test that compares the XML generated for 34 combinations of data
  and model choices with stored snapshots, and some 800 lines of commented-out code and
  seven unused classes were removed.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))
- The birth-death models moved to `dr.evomodel.birthdeath`; `EfficientSpeciationLikelihood`
  is now `EfficientBirthDeathLikelihood`.
- **MPI support has been removed**, along with the `-mpi` option, the `mpi.rank.postfix`
  property that labelled output files by rank, and the MPJ library. It had not worked for
  some time: the library was on the compile path but was never put into `beast.jar`, so a
  released BEAST could not have used it. Use `-threads`, the BEAGLE options, or `-mc3_chains`
  for parallel work.
- Removed `LogTransformedNormalKDEDistribution` and `LogitTransformedNormalKDEDistribution`,
  which duplicated about a thousand lines between them. `TransformedNormalKDEDistribution`
  has done the same job for any transform since 2019 and is what every caller already used;
  the two classes it replaced were left behind and did not give the right answer.
  ([#969](https://github.com/beast-dev/beast-mcmc/issues/969))
- `MultiPartitionTreeLikelihood` shares the count of BEAGLE instances with the other
  implementations, as they already share it with each other. BEAGLE devices are handed out
  by indexing `-beagle_order` with that count, so a class holding its own would start again
  at the first device. Nothing constructs this class, so this changes no behaviour; it keeps
  the bug fixed for v10.5.0 from returning with it.
  ([#1079](https://github.com/beast-dev/beast-mcmc/issues/1079))
- Eigen3 (JNI) and EJML eigendecomposition implementations were added but are not yet
  selectable from XML.

---

## BEAST X v10.5.0 *(2 July 2025)*

Initial public release of BEAST X (v10.5 series).

For earlier version history see the [BEAST website](https://beast.community/).
