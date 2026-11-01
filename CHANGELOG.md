# BEAST X Version History

Bayesian Evolutionary Analysis Sampling Trees  
https://beast.community/

All issues are tracked at https://github.com/beast-dev/beast-mcmc/issues

---

## BEAST X v10.6.0 *(in development)*

### New features

#### BEAST

- If the operator log file name has the extension `.csv` then a CSV file will be created (or if `format="csv"` is set).
- The operators element can now be logged to provide information about operator performance.
- **Output directory option** — `-output <directory>` on the command line, or "Write all
  output files to a folder" in the launch dialog, writes every log, tree, checkpoint and
  citation file to the chosen directory. The directory is checked to exist and be writable
  before the run starts. A checkpoint file named in the XML is now written next to the log
  files rather than in the working directory.
  ([#1267](https://github.com/beast-dev/beast-mcmc/pull/1267)) ([#1280](https://github.com/beast-dev/beast-mcmc/pull/1280))
- **Launch dialog** — the BEAGLE options are shown only when the BEAGLE library is found;
  otherwise the dialog explains that BEAGLE is required and where to get it. Quit and the
  window's close button now quit the application, and pressing Run without an XML file
  returns to the dialog instead of exiting.
- **Fast tree model** — `<treeModel implementation="fast">` keeps the tree in flat arrays
  stored copy-on-write. It produces the same MCMC output as the default tree model at a
  lower cost per step on large trees. Node rates and node traits are not supported by it.
  ([#1274](https://github.com/beast-dev/beast-mcmc/pull/1274))
- **Faster coalescent calculations** — the coalescent interval times are sorted
  incrementally rather than from scratch at every evaluation, and the CTMC scale prior keeps
  its tree length up to date incrementally.
  ([#1275](https://github.com/beast-dev/beast-mcmc/pull/1275))
- **Chain profiling** — `-Dmcmc.profile=steps`, `likelihoods` or `all` times the phases of
  every MCMC step per operator and every likelihood evaluation, reported after the operator
  analysis.
  ([#1276](https://github.com/beast-dev/beast-mcmc/pull/1276))
- **Integrated coalescent prior** — `IGCoalescentLikelihood` marginalises the constant
  population size under an inverse-gamma prior, ported from the work of Luiz Max Carvalho
  ([#1075](https://github.com/beast-dev/beast-mcmc/pull/1075)) and corrected for the shared population size.
  ([#1273](https://github.com/beast-dev/beast-mcmc/pull/1273))
- **Q.3Di.AF and Q.3Di.LLM substitution models** for FoldSeek's 3Di structural alphabet
  (Garg and Hochberg, 2025, Mol Biol Evol 42:msaf124), with their own equilibrium
  frequencies.
- **BASTA** — adjoint (reverse-mode) gradients for the structured coalescent migration
  rates, with BSSVS indicators applied through the chain rule, and support for
  travel-history (ancestral trait) trees.
- **Geodesic HMC** on spheres and Stiefel manifolds, with masks, and the correlation matrix
  of a factor model available as a parameter rather than only as a statistic.
- `duplicatedParameter` can take a fixed, constant number of copies as an attribute.
- `glmModel` takes an id that is used when logging its rates, and GLM coefficient and
  indicator dimensions are named after their predictors.
- On macOS, BEAST, TreeAnnotator and LogCombiner have the application menu (About, Quit)
  and the Dock icon, and BEAST's About box is laid out like BEAUti's.

#### BEAUti

- Added option to write doubling time in days to the log file when using a coalescent model with an exponential growth parameter.
- Added an option to write the epidemiological reproduction number (R0) to the log file when using a coalescent model with an exponential growth parameter. This requires an estimate of the mean and standard deviation of the generation time distribution to be specified in BEAUti.
- **Fixed local clock** — the Clocks panel shows the Yoder & Yang (2000) citation and lists
  the taxon sets that will be given their own rate (those marked monophyletic in the Taxa
  panel), or explains how to define one. The "Subtree prior model" drop-down is greyed out
  when no monophyletic taxon set exists.
- **Gao et al. (2023) prior specification for discrete trait models** — an option in the
  Sites panel replaces the default priors on a forward-in-time discrete trait model with the
  more diffuse specification of Gao, May, Rannala & Moore (2023, PNAS 120:e2213913120),
  including a new "Exponential with estimated mean" prior on the overall rate scalar.
- **macOS integration** — the application menu, the Dock icon and files opened from the
  Finder go through `java.awt.Desktop`, the title bar takes the window's colour, the toolbar
  icons are drawn as vectors for high-resolution displays, combo boxes use Aqua's pop-up
  style and highlight the item under the mouse, and `-Dswing.defaultlaf` overrides the
  system look and feel.
  ([#1277](https://github.com/beast-dev/beast-mcmc/pull/1277)) ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))
- Import errors are reported in a dialog that says what is wrong with the file, and the
  macOS file dialog accepts several files at once.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))

#### TreeAnnotator

- Can specify multiple HPD intervals e.g., using `--hpd 0.95,0.90,0.5'
- Can generate a KDE curve of the distribution of height values for visualisation in PearTree
- Clade MRCA heights can be calculated to avoid negative branch lengths (CA method of Heled and Bouckaert, (2013))

### Bug fixes

#### BEAST

- **Fixed a run resumed from a checkpoint not following the run that wrote it** — tips were
  numbered by their position in the (random) starting tree, so a resumed run had different
  node numbers and a permuted `branchRates.categories`, and adopting the saved tree permuted
  the categories again; relaxed clock analyses showed it first. Tips are now numbered by
  taxon name, in both the default and the fast tree model, and the saved tree is adopted
  without permuting node-indexed parameters. As a consequence a run started with the same
  seed as under an earlier version follows a different trajectory.
  ([#1225](https://github.com/beast-dev/beast-mcmc/issues/1225))

- **Fixed the branch-rate diagonal Hessian for transformed rates** — the reciprocal, exp
  and location-scale transforms gave a wrong mass matrix to HMC with (adaptive) diagonal
  preconditioning. Reporting a `traitDataLikelihood` with `arbitraryBranchRates` also no
  longer fails with an ArrayIndexOutOfBoundsException.
  ([#1279](https://github.com/beast-dev/beast-mcmc/pull/1279))

- **Fixed site pattern compression taking quadratic time** — a 70 taxon by 644,000 site
  alignment never finished loading in BEAST or BEAUti. Patterns are now found by a hash
  lookup.
  ([#1269](https://github.com/beast-dev/beast-mcmc/pull/1269))

- **Fixed BASTA giving -Inf for every state on an ancestral-trait tree**, where a
  travel-history ancestor sits on a zero-length coalescent interval by construction.

- **Fixed BigFastTreeIntervals** restoring stale intervals after a rejected proposal, and
  giving wrong intervals for tip-height moves after a full rebuild.

- **Fixed InverseGammaDistribution counting its normalising constant twice** in `pdf` and
  `logPdf`.

- **Fixed the latent liability model's variable sign** being taken from the starting value,
  which threw when a variable started exactly on zero; it is now taken from the bounds.

- **Fixed a NullPointerException in the rate covariance statistic** with the fast tree model.

#### BEAUti

- **Fixed BEAUti reverting clock-rate prior to CTMC reference prior** — custom
  rate priors (e.g. a lognormal prior on the substitution rate) specified in
  BEAUti were being silently replaced by the CTMC reference prior when the XML
  was written.
  ([#1238](https://github.com/beast-dev/beast-mcmc/issues/1238))

- Fixed an issue where if the name stem was too long, the names of the log files and tree files could be truncated in the XML.

- **Fixed saving and opening a BEAUti document** — saving failed for any analysis with an
  imported tree, and a document that did open had its data types replaced by copies so the
  Partitions and Site Models panels still showed the previous document.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))

- **Fixed the fixed local clock failing to generate XML** for a taxon set imported from a
  one-column file (a NullPointerException on the include-stem flag).
  ([#1231](https://github.com/beast-dev/beast-mcmc/issues/1231))

- **Fixed every taxon set being nested as a subtree of the random starting tree**, which
  made BEAST fail with "subtrees' taxa overlap" when sets overlapped. Only taxon sets marked
  monophyletic are nested now.

- **Fixed per-tree taxon lists not being written when partitions have different taxa**, so
  the XML failed with "Object with idref=X.taxa has not been previously declared".
  ([#1271](https://github.com/beast-dev/beast-mcmc/pull/1271))

- **Fixed the HMC relaxed clock and shrinkage local clock XML with several clock models.**
  ([#1272](https://github.com/beast-dev/beast-mcmc/pull/1272))

- **Fixed the gamma and Cauchy relaxed random walks** placing the gamma on a scalar of the
  diffusion variance rather than the precision, so the Cauchy setting was not a Cauchy.
  ([#1279](https://github.com/beast-dev/beast-mcmc/pull/1279))

- **Fixed a bit-flip operator on the rate indicators being written without BSSVS**, which
  referred to a parameter that was never written.
  ([#1278](https://github.com/beast-dev/beast-mcmc/pull/1278))

- **Fixed `includeStem` always written as true** for taxon set TMRCA statistics, whatever
  the checkbox said.

- **Fixed robust counting XML** for prefixed traits and for alignments with sites that are
  gaps in every sequence.

- **Fixed amino acid site model XML generation** failing with a NullPointerException.

- **Fixed importing a traits file turning on tip dates** for the whole analysis.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))

- **Fixed names not being escaped in the XML** — a partition, model, taxon set, trait or
  file name containing an ampersand, an angle bracket or a quote produced a file BEAST
  could not read.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))

- Generate no longer leaves a truncated XML file behind when it fails, takes the value
  being edited on the current panel before validating, and reports a bad combination of
  models before anything is written, switching to the panel concerned.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))

- Setting the Clock or Tree column of the Partitions table now follows the same rules as
  the link buttons, and unlinking or cloning a model copies the complete settings, so
  discrete and continuous trait models keep their structure.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))

- Dialogs are disposed when closed, imported files are closed, clearing a number field or
  selecting a separator in the tree prior list no longer throws, hovering over the MCMC
  panel's checkboxes no longer marks the document as changed, and declining to take a prior
  out of a hierarchical model abandons the edit.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))

- The stochastic Dollo model, which produced an XML whose likelihood took no part in the
  MCMC, has been removed.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))

- Taxon names containing a double quote are rejected on import, with the offending names
  listed, since such names cannot be written to a tree file.

#### TreeAnnotator

- **Fixed TreeAnnotator not summarising 2-dimensional locations correctly** —
  continuous phylogeographic location traits (e.g. `location1` / `location2`)
  were being renamed to `location.rate1` / `location.rate2` in the annotated
  tree output.
  ([#1242](https://github.com/beast-dev/beast-mcmc/issues/1242))

- **Fixed performance regression in TreeAnnotator "Collecting node information"
  stage** — processing of large tree sets became increasingly slow and eventually
  stalled, a regression introduced after Beta5.
  ([#1233](https://github.com/beast-dev/beast-mcmc/issues/1233))

- **Fixed bug in TreeAnnotator for continuous phylogeography with fixed trees** —
  coordinate annotations were being replaced by branch lengths and then discarded
  when the input tree topology was fixed rather than sampled.
  ([#1255](https://github.com/beast-dev/beast-mcmc/issues/1255))

- **Fixed a NullPointerException when a user target tree contains clades absent from the
  sampled trees** — such clades now get a posterior of 0 and keep their target-tree height.

### Other changes

- Updated HIPSTR citation to the final published paper.
  ([#1243](https://github.com/beast-dev/beast-mcmc/issues/1243))
- Taxon names that need quoting in tree and NEXUS files (whitespace or any of `. ( ) [ ] : ; , = '`)
  are now delimited with double quotes rather than single quotes in all of BEAST's output, so
  names containing an apostrophe can be written and read back. A taxon id containing a
  double quote is rejected when the XML is parsed.
- BEAGLE is loaded with `System.loadLibrary`, which JDK 24 and later warn about unless native
  access is enabled: the `beast.jar` manifest enables it for `java -jar`, and the Mac and
  Linux `beast` scripts add the option when the installed Java accepts it.
- Stopping a run uses `Thread.interrupt` rather than the removed `Thread.stop`.
- Development builds report their version as 10.x.x rather than an unexpanded `${version}`.
- BEAUti has a headless test that compares the XML generated for 34 combinations of data
  and model choices with stored snapshots, and some 800 lines of commented-out code and
  seven unused classes were removed.
  ([#1283](https://github.com/beast-dev/beast-mcmc/pull/1283))
- The Linux release includes the BEASTGen shell scripts and the Windows launcher configuration
  was updated.

---

## BEAST X v10.5.0 *(2 July 2025)*

Initial public release of BEAST X (v10.5 series).

For earlier version history see the [BEAST website](https://beast.community/).
