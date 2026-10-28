# Creates the brms models used by the tests (tests/testthat/fixtures/).
#
# Fitting brms models compiles Stan programs, which is too slow for the test
# suite, so small models are fitted once here and saved. The compiled Stan
# model is removed before saving; the draws, data and formula that
# relweights() and the tests use are kept.
#
# Run from the package root: Rscript data-raw/brms-fixtures.R

library(brms)

set.seed(20261028)
n <- 300
d <- data.frame(x1 = rnorm(n))
d$x2 <- 0.6 * d$x1 + rnorm(n, sd = 0.8)
d$x3 <- rnorm(n)
d$y <- 0.5 * d$x1 + 0.3 * d$x2 + 0.2 * d$x3 + rnorm(n)
d$y2 <- -0.4 * d$x1 + 0.5 * d$x2 + rnorm(n)

fit_brms <- function(formula) {
  brm(formula, data = d, chains = 2, iter = 600, warmup = 300, seed = 1,
      refresh = 0, save_warmup = FALSE, backend = "rstan")
}

strip <- function(fit) {
  fit$fit@stanmodel <- methods::new("stanmodel")
  fit$fit@.MISC <- new.env(parent = emptyenv())
  fit
}

dir.create("tests/testthat/fixtures", showWarnings = FALSE, recursive = TRUE)

uni <- fit_brms(bf(y ~ x1 + x2 + x3))
saveRDS(strip(uni), "tests/testthat/fixtures/brms-uni.rds", compress = "xz")

mv <- fit_brms(bf(mvbind(y, y2) ~ x1 + x2) + set_rescor(TRUE))
saveRDS(strip(mv), "tests/testthat/fixtures/brms-mv.rds", compress = "xz")
