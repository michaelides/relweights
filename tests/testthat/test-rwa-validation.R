# Validation against the 'rwa' package (Chan, CRAN), an independent
# implementation of Johnson's (2000) relative weights for data frames.

skip_if_not_installed("rwa")

rwa_multiple <- function(data, outcome, predictors, ...) {
  suppressMessages(rwa::rwa(data, outcome, predictors, method = "multiple",
                            sort = FALSE, ...))
}

expect_matches_rwa <- function(weights, r2, reference, tolerance = 1e-8) {
  result <- reference$result
  expect_equal(weights$predictor, as.character(result$Variables))
  expect_equal(weights$raw, result$Raw.RelWeight, tolerance = tolerance)
  expect_equal(weights$rescaled, result$Rescaled.RelWeight,
               tolerance = tolerance)
  expect_equal(r2, reference$rsquare, tolerance = tolerance)
  # rwa's Sign is the sign of the orthogonal-counterpart coefficient.
  expect_equal(ifelse(weights$sign_orth < 0, "-", "+"), result$Sign)
}

# Data whose correlation matrix equals the model-implied correlations exactly.
implied_data <- function(fit, variables, block = 1L, n = 500L) {
  R <- lavaan::lavInspect(fit, "cor.all")
  if (is.list(R)) R <- R[[block]]
  R <- unclass(R)[variables, variables]
  as.data.frame(MASS::mvrnorm(n, mu = rep(0, length(variables)), Sigma = R,
                              empirical = TRUE))
}

test_that("linear models match rwa", {
  rw <- relweights(lm(mpg ~ cyl + disp + hp + gear, data = mtcars))
  expect_matches_rwa(rw$weights, rw$r2$r2,
                     rwa_multiple(mtcars, "mpg", mtcars_predictors))

  predictors <- setdiff(names(mtcars), "mpg")
  rw <- relweights(lm(mpg ~ ., data = mtcars))
  expect_matches_rwa(rw$weights, rw$r2$r2,
                     rwa_multiple(mtcars, "mpg", predictors))
  # Here the two sign conventions disagree for disp.
  expect_equal(rw$weights$predictor[rw$weights$sign_beta != rw$weights$sign_orth],
               "disp")
})

test_that("weighted linear models match weighted rwa", {
  d <- mtcars
  d$w <- rep(c(0.5, 1, 1.5, 2.5), length.out = nrow(d))
  rw <- relweights(lm(mpg ~ cyl + disp + hp + gear, data = d, weights = w))
  expect_matches_rwa(rw$weights, rw$r2$r2,
                     rwa_multiple(d, "mpg", mtcars_predictors, weight = "w"))
})

test_that("observed-variable lavaan regressions match rwa", {
  fit <- quiet_sem("mpg ~ cyl + disp + hp + gear", data = mtcars)
  rw <- relweights(fit)
  expect_matches_rwa(rw$weights, rw$r2$r2,
                     rwa_multiple(mtcars, "mpg", mtcars_predictors),
                     tolerance = 1e-6)
})

test_that("latent regressions match rwa on their implied correlations", {
  skip_if_not_installed("MASS")

  fit <- quiet_sem(pd_model, data = lavaan::PoliticalDemocracy)
  rw <- relweights(fit, outcome = "dem65")
  d <- implied_data(fit, c("dem65", "ind60", "dem60"))
  expect_matches_rwa(rw$weights, rw$r2$r2,
                     rwa_multiple(d, "dem65", c("ind60", "dem60")))

  fit <- quiet_sem(hs_model, data = hs_data)
  rw <- relweights(fit)
  predictors <- c("visual", "textual", "ageyr", "grade")
  d <- implied_data(fit, c("speed", predictors))
  expect_matches_rwa(rw$weights, rw$r2$r2,
                     rwa_multiple(d, "speed", predictors))
})

test_that("multigroup regressions match rwa within each group", {
  fit <- quiet_sem("x9 ~ x7 + x8 + ageyr", data = hs_data, group = "school")
  rw <- relweights(fit)
  for (g in c("Pasteur", "Grant-White")) {
    reference <- rwa_multiple(hs_data[hs_data$school == g, ], "x9",
                              c("x7", "x8", "ageyr"))
    expect_matches_rwa(rw$weights[rw$weights$group == g, ],
                       rw$r2$r2[rw$r2$group == g], reference,
                       tolerance = 1e-6)
  }
})

test_that("multigroup latent regressions match rwa within each group", {
  skip_if_not_installed("MASS")

  fit <- quiet_sem(hs_model, data = hs_data, group = "school")
  rw <- relweights(fit)
  predictors <- c("visual", "textual", "ageyr", "grade")
  for (b in 1:2) {
    d <- implied_data(fit, c("speed", predictors), block = b)
    expect_matches_rwa(rw$weights[rw$weights$block == b, ],
                       rw$r2$r2[b], rwa_multiple(d, "speed", predictors))
  }
})

test_that("bootstrap replicates match rwa's replicates", {
  fit <- lm(mpg ~ cyl + disp + hp + gear, data = mtcars)
  set.seed(123)
  rw <- relweights(fit, bootstrap = TRUE, n_boot = 100)
  set.seed(123)
  reference <- suppressWarnings(
    rwa_multiple(mtcars, "mpg", mtcars_predictors, bootstrap = TRUE,
                 n_bootstrap = 100)
  )
  raw <- rw$boot$index$statistic == "raw"
  expect_equal(rw$boot$t[, raw], unname(reference$bootstrap$boot_object$t),
               tolerance = 1e-8)
  expect_equal(rw$boot$t0[raw], unname(reference$bootstrap$boot_object$t0),
               tolerance = 1e-8)
})
