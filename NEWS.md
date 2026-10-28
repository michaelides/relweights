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
