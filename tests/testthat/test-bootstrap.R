test_that("linear-model bootstrap adds percentile intervals", {
  fit <- lm(mpg ~ cyl + disp + hp + gear, data = mtcars)
  set.seed(1)
  rw <- relweights(fit, bootstrap = TRUE, n_boot = 200)
  w <- rw$weights

  expect_true(all(c("raw_lower", "raw_upper", "rescaled_lower",
                    "rescaled_upper") %in% names(w)))
  expect_true(all(w$raw_lower <= w$raw_upper))
  expect_true(all(w$rescaled_lower <= w$rescaled_upper))
  expect_true(all(w$raw_lower >= 0))
  expect_true(rw$r2$r2_lower <= rw$r2$r2_upper)
  expect_true(rw$r2$r2_lower < rw$r2$r2 && rw$r2$r2 < rw$r2$r2_upper)

  b <- rw$boot
  expect_equal(dim(b$t), c(200L, 9L))
  expect_equal(b$index$statistic, rep(c("raw", "rescaled", "r2"), c(4, 4, 1)))
  expect_equal(b$t0, c(w$raw, w$rescaled, rw$r2$r2))
  expect_equal(b$n_boot, 200L)
  expect_equal(b$n_failed, 0L)
  expect_equal(b$conf_level, 0.95)
  expect_identical(b$method, "boot::boot")
  expect_equal(rowSums(b$t[, 1:4]), b$t[, 9])
  expect_equal(unname(w$raw_lower[1]),
               unname(stats::quantile(b$t[, 1], 0.025, type = 6)))
})

test_that("bootstrap results are reproducible with set.seed()", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  set.seed(42)
  first <- relweights(fit, bootstrap = TRUE, n_boot = 50)
  set.seed(42)
  second <- relweights(fit, bootstrap = TRUE, n_boot = 50)
  expect_identical(first$weights, second$weights)
  expect_identical(first$boot$t, second$boot$t)
})

test_that("the confidence level changes the intervals", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  set.seed(42)
  wide <- relweights(fit, bootstrap = TRUE, n_boot = 100, conf_level = 0.99)
  set.seed(42)
  narrow <- relweights(fit, bootstrap = TRUE, n_boot = 100, conf_level = 0.5)
  expect_true(all(wide$weights$raw_lower <= narrow$weights$raw_lower))
  expect_true(all(wide$weights$raw_upper >= narrow$weights$raw_upper))
})

test_that("weighted models resample rows together with their weights", {
  d <- mtcars
  d$w <- rep(c(1, 2, 4), length.out = nrow(d))
  fit <- lm(mpg ~ wt + hp, data = d, weights = w)
  set.seed(3)
  rw <- relweights(fit, bootstrap = TRUE, n_boot = 30)
  # The R2 of each replicate equals the R2 of the weighted regression refitted
  # to the same resampled rows.
  set.seed(3)
  reference <- boot::boot(d, function(data, i) {
    summary(lm(mpg ~ wt + hp, data = data[i, ], weights = w))$r.squared
  }, R = 30)
  expect_equal(rw$boot$t[, 5], reference$t[, 1])
})

test_that("replicates that cannot be analysed are excluded with a warning", {
  d <- mtcars
  d$rare <- c(1, rep(0, nrow(d) - 1L))
  fit <- lm(mpg ~ wt + rare, data = d)
  set.seed(1)
  expect_warning(rw <- relweights(fit, bootstrap = TRUE, n_boot = 60),
                 "could not be analysed")
  expect_gt(rw$boot$n_failed, 0L)
  failed <- !stats::complete.cases(rw$boot$t)
  expect_equal(sum(failed), rw$boot$n_failed)
  expect_equal(unname(rw$weights$raw_upper[1]),
               unname(rw_percentile(rw$boot$t[!failed, 1], 0.95)[2]))
  expect_output(print(rw), "failed\\s+and\\s+were\\s+excluded")
})

test_that("bootstrap arguments are checked", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_error(relweights(fit, bootstrap = TRUE, n_boot = 0), "n_boot")
  expect_error(relweights(fit, bootstrap = TRUE, n_boot = 10.5), "n_boot")
  expect_error(relweights(fit, bootstrap = TRUE, n_boot = "10"), "n_boot")
  expect_error(relweights(fit, bootstrap = TRUE, conf_level = 1), "conf_level")
  expect_error(relweights(fit, bootstrap = TRUE, conf_level = c(0.9, 0.95)),
               "conf_level")
  expect_error(relweights(fit, bootstrap = NA), "TRUE or FALSE")
  expect_error(relweights(fit, bootstrap = TRUE, iseed = 1),
               "Unknown argument\\(s\\).*'iseed'")
  set.seed(1)
  expect_no_error(relweights(fit, bootstrap = TRUE, n_boot = 10,
                             parallel = "no", ncpus = 1))
})

test_that("lavaan models are bootstrapped by refitting", {
  fit <- quiet_sem("x9 ~ x7 + x8 + ageyr", data = hs_data)
  set.seed(10)
  rw <- relweights(fit, bootstrap = TRUE, n_boot = 20)

  expect_equal(dim(rw$boot$t), c(20L, 7L))
  expect_identical(rw$boot$method, "lavaan::lavBootstrap")
  expect_equal(rw$boot$n_failed, 0L)
  expect_true(all(rw$weights$raw_lower <= rw$weights$raw_upper))
  expect_true(rw$r2$r2_lower <= rw$r2$r2_upper)
  expect_equal(rw$boot$t0, c(rw$weights$raw, rw$weights$rescaled, rw$r2$r2))

  set.seed(10)
  again <- relweights(fit, bootstrap = TRUE, n_boot = 20)
  expect_identical(again$boot$t, rw$boot$t)
})

test_that("multigroup lavaan models are bootstrapped within groups", {
  fit <- quiet_sem("x9 ~ x7 + x8", data = hs_data, group = "school")
  set.seed(11)
  rw <- relweights(fit, bootstrap = TRUE, n_boot = 20, type = "ordinary")
  expect_equal(dim(rw$boot$t), c(20L, 10L))
  expect_false(anyNA(rw$r2$r2_lower))
  expect_equal(rw$boot$index$regression, rep(1:2, each = 5))
})

test_that("latent-variable models can be bootstrapped", {
  fit <- quiet_sem(hs_model, data = hs_data)
  set.seed(12)
  rw <- suppressWarnings(relweights(fit, bootstrap = TRUE, n_boot = 20))
  expect_equal(nrow(rw$boot$t), 20L)
  expect_lt(rw$boot$n_failed, 20L)
  expect_true(all(rw$weights$raw_lower <= rw$weights$raw_upper))
})

test_that("unsupported lavaan bootstraps are informative errors", {
  multilevel <- quiet_sem(twolevel_model, data = lavaan::Demo.twolevel,
                          cluster = "cluster")
  expect_error(relweights(multilevel, bootstrap = TRUE, n_boot = 5),
               "not available for multilevel")

  S <- cov(mtcars[, c("mpg", mtcars_predictors)])
  moments <- quiet_sem("mpg ~ cyl + disp + hp + gear", sample.cov = S,
                       sample.nobs = 32)
  expect_error(relweights(moments, bootstrap = TRUE, n_boot = 5), "raw data")
})
