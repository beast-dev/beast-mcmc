# BEAST X Version History

Bayesian Evolutionary Analysis Sampling Trees  
https://beast.community/

All issues are tracked at https://github.com/beast-dev/beast-mcmc/issues

---

## BEAST X v10.5.1 *(in development)*

### Bug fixes

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

#### BEAUti

- **Fixed BEAUti reverting clock-rate prior to CTMC reference prior** — custom
  rate priors (e.g. a lognormal prior on the substitution rate) specified in
  BEAUti were being silently replaced by the CTMC reference prior when the XML
  was written.
  ([#1238](https://github.com/beast-dev/beast-mcmc/issues/1238))

### Other changes

- Updated HIPSTR citation to the final published paper.
  ([#1243](https://github.com/beast-dev/beast-mcmc/issues/1243))

---

## BEAST X v10.5.0 *(2 July 2025)*

Initial public release of BEAST X (v10.5 series).

For earlier version history see the [BEAST website](https://beast.community/).
