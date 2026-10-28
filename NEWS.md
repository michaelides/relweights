# relweights 0.2.0

* `relweights()` now analyses Bayesian models, computing the weights in every
  posterior draw:
  * linear models fitted with brms (gaussian family, identity link,
    population-level effects only; multivariate models one response at a
    time). R-squared is the model-based Bayesian R-squared of Gelman et al.
    (2019); `r2_model` holds the residual-based value of `brms::bayes_R2()`.
  * structural equation models fitted with blavaan, including multigroup
    and two-level models. The model-implied correlations are recomputed in
    every draw.
* Results of Bayesian models report posterior means with equal-tailed
  credible intervals, the posterior probabilities that `beta` and `orth` are
  positive, and R-hat and bulk ESS of the raw weights. The per-draw results
  are kept in `$draws`, and the new `ndraws` argument thins the draws.
* `as.data.frame()` gains `what = "ranks"` and `what = "pairwise"`, which
  return posterior rank probabilities and the probabilities that one weight
  exceeds another.
* `print()` and the plots label credible intervals as such.
* New vignette: `vignette("bayesian")`.

# relweights 0.1.0

* First release.
* `relweights()` computes Johnson's (2000) relative weights for linear models
  fitted with `lm()`, including models with prior weights. It also analyses
  every structural regression of a model fitted with lavaan, using the
  model-implied correlations among observed and latent variables, separately
  for each group and level.
* Results include raw and rescaled weights, standardized coefficients, the
  coefficients of the orthogonal counterparts, zero-order correlations and
  two sign conventions (`sign_beta` and `sign_orth`).
* `bootstrap = TRUE` adds percentile intervals for the raw weights, the
  rescaled weights and R-squared. Linear models are resampled with
  `boot::boot()`, and lavaan models are refitted with
  `lavaan::lavBootstrap()`.
* `autoplot()` and `plot()` draw four ggplot2 plots: the weights, the
  decomposition of R-squared, a comparison of importance indices, and a
  comparison across groups and levels.
* The test suite checks the results against the rwa package.
