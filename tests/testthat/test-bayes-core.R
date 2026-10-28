lm_fit <- lm(mpg ~ cyl + disp + hp + wt, data = mtcars)
lm_X <- model.matrix(lm_fit)[, -1]
lm_y <- mtcars$mpg

# Synthetic "posterior draws" around the least-squares estimates.
synthetic_draws <- function(n = 60) {
  set.seed(3)
  b <- coef(lm_fit)[-1]
  B <- matrix(b, n, length(b), byrow = TRUE) *
    matrix(stats::runif(n * length(b), 0.5, 1.5), n)
  sigma <- stats::runif(n, 2, 3)
  rw_fixed_x_draws(B, sigma, lm_X, lm_y, "mpg")
}

synthetic_result <- function(draws = synthetic_draws(), chains = 2L) {
  n <- length(draws$r2)
  spec <- list(block = 1L, group = NA_character_, level = NA_character_,
               outcome = "mpg", predictors = colnames(lm_X),
               terms = colnames(lm_X), n = nrow(lm_X), r2_model = NA_real_,
               label = "mpg")
  info <- list(model_class = "brmsfit", call = NULL, weighted = FALSE,
               ngroups = 1L, nlevels = 1L, group_labels = character(),
               level_labels = character(), multi_column_terms = character(),
               formula = "mpg ~ cyl + disp + hp + wt",
               r2_model_label = "bayes_R2")
  chain <- rep(seq_len(chains), each = n / chains)
  iteration <- rep(seq_len(n / chains), chains)
  rw_bayes_build(list(spec), list(draws), chain, iteration, 0.9, info)
}

test_that("least-squares estimates fed to the per-draw code reproduce lm()", {
  sigma <- sqrt(sum(residuals(lm_fit)^2) / (nrow(lm_X) - 1))
  out <- rw_fixed_x_draws(matrix(coef(lm_fit)[-1], nrow = 1), sigma, lm_X,
                          lm_y, "mpg")
  ref <- relweights(lm_fit)

  for (column in c("raw", "rescaled", "beta", "orth", "r")) {
    expect_equal(as.vector(out[[column]]), ref$weights[[column]],
                 tolerance = 1e-10)
  }
  expect_equal(out$r2, summary(lm_fit)$r.squared, tolerance = 1e-10)
  expect_equal(out$r2_model, summary(lm_fit)$r.squared, tolerance = 1e-10)
  expect_equal(out$lambda %*% out$lambda, cor(lm_X), ignore_attr = TRUE)
})

test_that("every draw matches the analysis of its own correlation matrix", {
  out <- synthetic_draws()
  vars <- c("mpg", colnames(lm_X))
  for (d in c(1, 30, 60)) {
    r <- out$r[d, ]
    R <- rbind(c(1, r), cbind(r, cor(lm_X)))
    dimnames(R) <- list(vars, vars)
    ref <- rw_compute(R, "mpg", colnames(lm_X))
    expect_equal(out$raw[d, ], ref$raw)
    expect_equal(out$rescaled[d, ], ref$rescaled)
    expect_equal(out$beta[d, ], ref$beta)
    expect_equal(out$orth[d, ], ref$orth)
    expect_equal(out$r2[d], ref$r2)
  }
})

test_that("posterior means are additive and intervals are quantiles", {
  draws <- synthetic_draws()
  x <- synthetic_result(draws)
  w <- x$weights

  expect_equal(w$raw, unname(colMeans(draws$raw)))
  expect_equal(w$beta, unname(colMeans(draws$beta)))
  expect_equal(sum(w$raw), x$r2$r2)
  expect_equal(sum(w$rescaled), 100)
  expect_equal(x$r2$r2, mean(draws$r2))
  expect_equal(x$r2$r2_model, mean(draws$r2_model))
  expect_equal(w$raw_lower[2],
               unname(stats::quantile(draws$raw[, 2], 0.05)))
  expect_equal(x$r2$r2_upper, unname(stats::quantile(draws$r2, 0.95)))
  expect_equal(w$p_beta_pos, unname(colMeans(draws$beta > 0)))
  expect_equal(w$p_orth_pos, unname(colMeans(draws$orth > 0)))
  expect_true(all(is.finite(w$rhat)))
  expect_true(all(w$ess_bulk > 0))
  expect_identical(x$info$intervals, "credible")
  expect_equal(x$info$ndraws, 60L)
  expect_equal(x$info$nchains, 2L)
  expect_named(x$draws$regressions, "mpg")
  expect_equal(x$draws$regressions$mpg$raw, draws$raw)
  expect_null(x$boot)
})

test_that("rank and pairwise probabilities are consistent", {
  draws <- synthetic_draws()
  x <- synthetic_result(draws)
  ranks <- as.data.frame(x, what = "ranks")
  p <- ncol(lm_X)

  expect_equal(nrow(ranks), p * p)
  expect_equal(as.vector(tapply(ranks$probability, ranks$predictor, sum)),
               rep(1, p))
  expect_equal(as.vector(tapply(ranks$probability, ranks$rank, sum)),
               rep(1, p))
  best <- colnames(lm_X)[max.col(draws$raw, ties.method = "first")]
  first <- ranks[ranks$rank == 1, ]
  expect_equal(first$probability,
               as.vector(table(factor(best, levels = colnames(lm_X)))) / 60)
  expect_equal(rw_first_rank(x), first$probability)

  pairs <- as.data.frame(x, what = "pairwise")
  expect_equal(nrow(pairs), p * (p - 1))
  forward <- pairs$probability[pairs$predictor == "wt" & pairs$other == "hp"]
  backward <- pairs$probability[pairs$predictor == "hp" & pairs$other == "wt"]
  expect_equal(forward, mean(draws$raw[, "wt"] > draws$raw[, "hp"]))
  expect_equal(forward + backward, 1)
})

test_that("failed draws are excluded with a warning", {
  draws <- synthetic_draws()
  draws$raw[c(2, 5), ] <- NA
  draws$r2[c(2, 5)] <- NA
  expect_warning(x <- synthetic_result(draws), "2 of 60 posterior draws")
  expect_equal(x$draws$n_failed, 2L)
  expect_equal(x$info$ndraws, 58L)
  expect_equal(x$weights$raw, unname(colMeans(draws$raw[-c(2, 5), ])))
  expect_true(all(is.na(x$weights$rhat)))

  draws$raw[] <- NA
  expect_error(suppressWarnings(synthetic_result(draws)),
               "None of the 60 posterior draws")
})

test_that("draws are thinned evenly within chains", {
  expect_equal(rw_thin_iterations(100, 2, NULL), 1:100)
  expect_equal(rw_thin_iterations(100, 2, 500), 1:100)
  expect_equal(rw_thin_iterations(100, 2, 10), c(1, 26, 50, 75, 100))
  expect_equal(rw_thin_iterations(100, 4, 2), 1)
})

test_that("Bayesian arguments are checked", {
  expect_error(rw_check_bayes_args(NULL, 0.95, list(bootstrap = TRUE)),
               "Bootstrap intervals are not used")
  expect_error(rw_check_bayes_args(NULL, 0.95, list(n_boot = 10)),
               "Bootstrap intervals are not used")
  expect_error(rw_check_bayes_args(NULL, 0.95, list(foo = 1)),
               "Unknown argument\\(s\\).*'foo'")
  expect_error(rw_check_bayes_args(NULL, 0.95, list(1)), "'<unnamed>'")
  expect_error(rw_check_bayes_args(0, 0.95, list()), "ndraws")
  expect_error(rw_check_bayes_args(2.5, 0.95, list()), "ndraws")
  expect_error(rw_check_bayes_args(NULL, 1, list()), "conf_level")
  expect_true(rw_check_bayes_args(10, 0.9, list()))
})

test_that("rank and pairwise tables need posterior draws", {
  rw <- relweights(lm_fit)
  expect_error(as.data.frame(rw, what = "ranks"), "posterior draws")
  expect_error(as.data.frame(rw, what = "pairwise"), "posterior draws")
})
