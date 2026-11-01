# BEAST X v10.6.0 Release Notes

BEAST X v10.6.0 adds new substitution and coalescent models, faster MCMC on
large trees, an output directory option, and many improvements to BEAUti, including the fixed local clock and new
prior specifications for discrete trait models. A full list of changes and
bug fixes is in [CHANGELOG.md](CHANGELOG.md).

## New models

- **Q.3Di substitution models for protein structure.** The `Q.3Di.AF` and
  `Q.3Di.LLM` amino acid models describe the evolution of FoldSeek's 3Di
  structural alphabet, in which each residue is coded by its local tertiary
  interactions rather than its identity (Garg and Hochberg, 2025, *Mol Biol
  Evol* 42:msaf124). Each comes with its own equilibrium frequencies and is
  selected like any other empirical amino acid model.

- **Integrated coalescent prior.** The `coalescentLikelihoodIG` element gives
  the constant-size coalescent with the population size integrated out under
  an inverse-gamma prior, so the population size does not have to be sampled.
  This follows the work of Luiz Max Carvalho and is useful as a tree prior
  for many short trees or where the population size itself is not of
  interest.

- **BASTA structured coalescent improvements.** BASTA can now provide
  adjoint (reverse-mode) gradients of the migration rates, including through
  BSSVS indicator variables, which allows the rates to be sampled with
  Hamiltonian Monte Carlo. It also supports trees with travel-history
  (ancestral trait) nodes.

- **Geodesic Hamiltonian Monte Carlo.** The
  `geodesicHamiltonianMonteCarloOperator` samples parameters that live on a
  sphere or a Stiefel manifold, such as orthonormal loadings, with optional
  masks. The correlation matrix of a factor model is now available as a
  parameter so it can be given a prior and sampled directly.

- **Gao et al. (2023) prior specification for discrete trait models.** In
  BEAUti's Sites panel the default priors of a forward-in-time discrete trait
  model can be replaced with the more diffuse specification of Gao, May,
  Rannala and Moore (2023, *PNAS* 120:e2213913120), including a new
  "Exponential with estimated mean" prior on the overall rate scalar.

## New features in BEAST

- **Output directory.** All output files (logs, trees, checkpoints and the
  citations file) can be written to a chosen directory with the
  `-output <directory>` command line option, or with the "Write all output
  files to a folder" option in the launch dialog. The directory is checked
  before the run starts. A checkpoint file named in the XML is written next
  to the log files rather than in the working directory.

- **Speed improvements.** A new tree model implementation
  (`<treeModel implementation="fast">`) and incremental coalescent
  calculations make each MCMC step cheaper on large trees, with identical
  results to earlier versions.

- **Operator logging.** The operators element can be logged to provide
  information about operator performance during a run. If the operator log
  file name ends in `.csv` (or `format="csv"` is given) the log is written as
  a CSV file.

- **Chain profiling.** Running with `-Dmcmc.profile=steps`, `likelihoods`
  or `all` times the phases of every MCMC step, per operator, and every
  likelihood evaluation. The timings are reported after the operator analysis
  at the end of the run.

- **Reproducible checkpointing.** A run resumed from a checkpoint now follows
  the run that wrote it exactly, including for relaxed clock models
  ([#1225](https://github.com/beast-dev/beast-mcmc/issues/1225)). Tips are
  numbered by taxon name rather than by their position in the random
  starting tree, so note that a run started with the same seed as under an
  earlier version will follow a different trajectory.

- **Launch dialog.** The BEAGLE options are shown only when the BEAGLE
  library is found; otherwise the dialog explains that BEAGLE is required and
  where to get it. Quit and the window's close button quit the application,
  and pressing Run without an XML file returns to the dialog.

- **XML conveniences.** `duplicatedParameter` can take a fixed number of
  copies as an attribute, `glmModel` takes an id that is used when logging
  its rates, and GLM coefficient and indicator dimensions are named after
  their predictors in the log.

- **macOS integration.** BEAST, TreeAnnotator and LogCombiner have the
  standard application menu (About, Quit) and Dock icon.

## New features in BEAUti

- **Fixed local clock.** The Clocks panel supports the fixed local clock of
  Yoder and Yang (2000), giving a separate rate to each taxon set marked
  monophyletic in the Taxa panel, and shows which taxon sets will receive
  their own rate.

- **Epidemiological quantities from exponential growth.** When a coalescent
  model with an exponential growth rate is used, the doubling time in days
  and the reproduction number R0 can be written to the log file. R0 requires
  the mean and standard deviation of the generation time distribution to be
  specified in BEAUti.

- **Gao et al. (2023) prior specification** for discrete trait models, as
  described under new models above.

- **Import error reporting.** Errors in imported files are reported in a
  dialog that says what is wrong, including taxon names containing a double
  quote, which cannot be written to a tree file and are now rejected with the
  offending names listed.

- **Safer XML generation.** Generate validates the value being edited on the
  current panel, reports an incompatible combination of models before
  anything is written and switches to the panel concerned, and never leaves a
  truncated XML file behind. Names containing an ampersand, angle bracket or
  quote are escaped correctly.

- **macOS integration.** The application menu, Dock icon and files opened
  from the Finder are handled natively, the title bar takes the window's
  colour, toolbar icons are drawn as vectors for high-resolution displays,
  combo boxes use the Aqua pop-up style, and the macOS file dialog accepts
  several files at once. The `-Dswing.defaultlaf` system property overrides
  the system look and feel.

- **Saving and opening documents** works again for analyses with imported
  trees, and the stochastic Dollo model, whose likelihood took no part in the
  MCMC, has been removed.

## New features in TreeAnnotator

- **Node height distributions for PearTree.** TreeAnnotator can write a
  kernel density estimate (KDE) of the posterior distribution of each node's
  height into the annotated tree, so that the full distribution, rather than
  only an HPD interval, can be visualised in PearTree.

- **Multiple HPD intervals.** Several HPD levels can be annotated at once,
  for example `--hpd 0.95,0.90,0.5`.

- **Common ancestor heights.** Clade MRCA heights can be calculated with the
  CA method of Heled and Bouckaert (2013) to avoid negative branch lengths.

- **Target trees with unsampled clades.** A user-supplied target tree
  containing clades absent from the sampled trees is annotated with a
  posterior of 0 for those clades, keeping the target-tree height, rather
  than failing.

## Other improvements

- **Double-quoted taxon names.** Taxon names that need quoting in tree and
  NEXUS files (whitespace or any of `. ( ) [ ] : ; , = '`) are delimited with
  double quotes rather than single quotes in all of BEAST's output, so names
  containing an apostrophe can be written and read back. A taxon id
  containing a double quote is rejected when the XML is parsed.

- **Large alignments load quickly.** Site pattern compression took time
  quadratic in the number of sites, so a 70 taxon by 644,000 site alignment
  never finished loading. Patterns are now found by a hash lookup.

- **Java 24 and later.** The warning about native access when BEAGLE is
  loaded is suppressed, by the `beast.jar` manifest when running with
  `java -jar` and by the Mac and Linux launch scripts.

- **Platform changes.** Stopping a run no longer relies on the removed
  `Thread.stop`, the Linux release includes the BEASTGen shell scripts, and
  the Windows launcher configuration has been updated.

- **Citations.** The HIPSTR citation has been updated to the final published
  paper.

Many bug fixes are included in this release, notably to the branch-rate
Hessian used by HMC relaxed clocks, the gamma and Cauchy relaxed random walks,
BASTA on travel-history trees, and the summarising of continuous
phylogeographic locations in TreeAnnotator. See
[CHANGELOG.md](CHANGELOG.md) for the complete list.
