# Site rate categories in the MAT likelihood

A discrete site rate model (e.g. `<siteModel>` with `<gammaShape gammaCategories="4">`) is supported by
the MAT likelihood with every branch history model except a per-site rate parameter (`siteRates`), and by
the history resampler's local and path blocks. Implemented in 19fe87a2e (likelihood, store, category
operator) and 8fed6a989 (resampler).

## Model

BEAST's site rate models are mixtures: a site's likelihood is the sum over categories k of p_k L(site | r_k).
Summing the categories out of a site with an explicit history would tie every site's term to every
branch, so the MAT augments instead: each variable site's category c_p is part of the sampled state, as
its history is. Given the categories, the likelihood decomposes by branch exactly as before, each site at
rate r_{c_p}, and the prior of the categories adds sum_p log p_{c_p}. Summed over the categories, the
posterior is that of the site rate model.

Constant patterns are pooled by state (constantOnly compression), so they cannot carry categories. They
have no events, so they are summed over the categories analytically: in category k and state a a constant
site's likelihood is exp of the sum over branches of the history model's unchanged-site term
u_a(r_k t_b), so the pooled patterns contribute W_a log sum_k p_k exp(sum_b u_a(r_k t_b)). The delegate
keeps the per-branch terms and their running totals, updated as branch terms are.

## Pieces

- **Store** (`DefaultMutationStore`): each pattern's category (journaled), and each node's category
  composition, the summed weights of its variable patterns in each (category, state), kept up to date by
  `setState` and `setCategory`, and discarded for the nodes whose states changed on revert (recomputed on
  demand).
- **Delegate**: pattern rates from the categories (constant patterns get rate 0, so exit-rate sums cover
  only the variable sites); `getCategoryBranchLogLikelihood` on the history models (integrated hits uses the
  category composition; the exit-rate models need nothing more); the constant-site and category-prior
  terms.
- **`siteRateCategoryGibbsOperator`**: draws every variable site's category from p_k L(history_p | r_k).
  A site's log likelihood at rate r is that of its root state staying unchanged on every branch,
  D_s(root; r), plus at each of its events on the branch above node n (from a to b) the event's term minus
  the unchanged term, minus D_a(n; r) plus D_b(n; r), where D_s(n; r) is the sum of the unchanged-site
  terms in state s over the branches below n (an event changes the state of the whole subtree below it,
  until the next event). A sweep costs O(nodes x categories x states + events x categories), plus a
  full recomputation of the likelihood. The categories mix quickly, so a weight of about 1% of the subtree
  moves' weight is enough: on small-2, 10% (16.5) took 1.0 minute per million states against 0.56 at 1%.
- **Resampler**: the categories are its rate classes. The sparse passes count the uniform sites, and
  group the ambiguity-only sites, by (category, state), each with its own conditional, so moves still cost
  O(events + changes). The block keeps branch densities and linear matrices per rate class.

## Cost

With integrated hits the constant sites' and the branch terms' unchanged parts need only the diagonals
log P_aa(r_c t), computed from the rate matrix's eigen decomposition with BaseSubstitutionModel's arithmetic
(the same values as the transition-probability cache); they are versioned by the rates, so a category sweep
does not recompute them. A sweep gives the store all the categories and the nodes' category compositions at
once (each node's from its parent's and the events on its branch), and the compositions are journaled. With
4 categories: small-2 0.33 minutes per million states (0.18 without), big-3 1.56 (0.76 without).

## Validation

- Unit tests: the likelihood against a direct site-by-site computation (through category, height and
  topology changes and restores), the category draws against the brute-force conditional, the local
  block's Hastings ratio against enumeration, and the local block's and the in-place Gibbs update's draws
  against the exact conditional, with single and integrated hits; each fails against a wrong reference.
- MCMC full-evaluation check at every state on small-2 with gamma categories.
- small-2 with a discrete gamma (4 categories) against BEAGLE:
  `optimization_notes/mat_local_spr/benchmark/make_validation_xmls.py --gamma` (see the benchmark README).
