skip_if_not_installed("brms")
skip_if_not_installed("posterior")

test_that("per-draw R2 values match brms", {
  fit <- brms_fixture("uni")
  rw <- relweights(fit)
  draws <- rw$draws$regressions$y

  # The residual-based R2 reproduces brms::bayes_R2() draw by draw.
  expect_equal(draws$r2_model, unname(brms::bayes_R2(fit, summary = FALSE)[, 1]),
               tolerance = 1e-10)
  # The model-based R2 is var(fitted) / (var(fitted) + sigma^2).
  epred <- brms::posterior_epred(fit)
  sigma <- as.vector(posterior::as_draws_matrix(fit, variable = "sigma"))
  var_fit <- apply(epred, 1, stats::var)
  expect_equal(draws$r2, var_fit / (var_fit + sigma^2), tolerance = 1e-10)
  expect_equal(rowSums(draws$raw), draws$r2)
  expect_equal(rw$r2$r2, mean(draws$r2))
  expect_equal(rw$r2$r2_model, mean(draws$r2_model))
})

test_that("posterior means are close to the least-squares results", {
  fit <- brms_fixture("uni")
  rw <- relweights(fit)
  ref <- relweights(lm(y ~ x1 + x2 + x3, data = fit$data))

  expect_equal(rw$weights$predictor, c("x1", "x2", "x3"))
  expect_equal(rw$r2$n, 300L)
  expect_lt(max(abs(rw$weights$rescaled - ref$weights$rescaled)), 5)
  expect_lt(abs(rw$r2$r2 - ref$r2$r2), 0.03)
  expect_lt(max(abs(rw$weights$beta - ref$weights$beta)), 0.05)
  expect_equal(sum(rw$weights$rescaled), 100)
  expect_true(all(rw$weights$raw_lower < rw$weights$raw &
                    rw$weights$raw < rw$weights$raw_upper))
  expect_true(all(rw$weights$p_beta_pos > 0.99))
  expect_true(all(rw$weights$rhat < 1.05))
  expect_identical(rw$info$model_class, "brmsfit")
  expect_identical(rw$info$formula, "y ~ x1 + x2 + x3")
})

test_that("multivariate models are analysed one response at a time", {
  fit <- brms_fixture("mv")
  rw <- relweights(fit)
  r2 <- brms::bayes_R2(fit, summary = FALSE)

  expect_equal(rw$r2$outcome, c("y", "y2"))
  expect_equal(rw$weights$outcome, rep(c("y", "y2"), each = 2))
  expect_equal(rw$draws$regressions$y$r2_model, unname(r2[, "R2y"]),
               tolerance = 1e-10)
  expect_equal(rw$draws$regressions$y2$r2_model, unname(r2[, "R2y2"]),
               tolerance = 1e-10)
  # x1 has a negative coefficient for y2.
  expect_lt(rw$weights$p_beta_pos[rw$weights$outcome == "y2"][1], 0.01)

  only <- relweights(fit, outcome = "y2")
  expect_equal(only$r2$outcome, "y2")
  expect_equal(only$weights$raw, rw$weights$raw[3:4])
  expect_error(relweights(fit, outcome = "y3"),
               "'y3' are not responses.*'y', 'y2'")
})

test_that("ndraws thins the draws within each chain", {
  fit <- brms_fixture("uni")
  rw <- relweights(fit, ndraws = 100)
  expect_equal(rw$info$ndraws, 100L)
  expect_equal(rw$draws$chain, rep(1:2, each = 50))
  expect_equal(rw$draws$iteration[1:3], c(1, 7, 13))

  one <- suppressWarnings(relweights(fit, ndraws = 1))
  expect_equal(one$info$ndraws, 2L)
  expect_true(all(is.na(one$weights$rhat)))
})

test_that("print, tables and plots describe credible intervals", {
  fit <- brms_fixture("uni")
  rw <- relweights(fit, conf_level = 0.9)
  out <- paste(capture.output(print(rw)), collapse = " ")
  out <- gsub("\\s+", " ", out)

  expect_match(out, "Model: brms model (gaussian, identity link): y ~ x1 + x2 + x3",
               fixed = TRUE)
  expect_match(out, "raw CrI", fixed = TRUE)
  expect_match(out, "(bayes_R2: ", fixed = TRUE)
  expect_match(out, "P(rank 1)", fixed = TRUE)
  expect_match(out, "Estimates are posterior means", fixed = TRUE)
  expect_match(out, "90% equal-tailed credible intervals from 600 posterior draws (2 chains)",
               fixed = TRUE)
  expect_match(out, "largest R-hat", fixed = TRUE)

  ranks <- as.data.frame(rw, what = "ranks")
  expect_equal(nrow(ranks), 9L)
  pairs <- as.data.frame(rw, what = "pairwise")
  expect_equal(nrow(pairs), 6L)

  for (type in c("weights", "r2", "compare")) {
    p <- autoplot(rw, type = type, sign = "beta")
    expect_no_warning(ggplot2::ggplot_build(p))
  }
  p <- autoplot(rw)
  expect_true("GeomErrorbar" %in% layer_geoms(p))
  expect_match(p$labels$caption, "90% equal-tailed credible intervals")
  expect_match(p$labels$caption, "posterior means")
  expect_false("GeomErrorbar" %in% layer_geoms(autoplot(rw, ci = FALSE)))
  expect_error(autoplot(rw, type = "groups"), "single group and level")
})

test_that("unsupported brms models are informative errors", {
  expect_error(relweights(empty_brm(y ~ x + (1 | g))), "group-level terms")
  expect_error(relweights(empty_brm(brms::bf(y ~ x, sigma ~ z))),
               "distributional parameter\\(s\\) 'sigma'")
  expect_error(relweights(empty_brm(y ~ s(x))), "smooth terms")
  expect_error(relweights(empty_brm(y ~ mo(o))), "special terms")
  expect_error(relweights(empty_brm(y ~ x + offset(z))), "an offset")
  expect_error(relweights(empty_brm(y | weights(w) ~ x)), "addition terms")
  expect_error(
    relweights(empty_brm(brms::bf(y ~ a * exp(x), a ~ 1, nl = TRUE))),
    "non-linear"
  )
  expect_error(relweights(empty_brm(b ~ x, family = brms::bernoulli())),
               "gaussian family and identity link")
  expect_error(relweights(empty_brm(y ~ 0 + x)), "no intercept")
  expect_error(relweights(empty_brm(y ~ x)), "no posterior draws")
  expect_warning(
    expect_error(relweights(empty_brm(y ~ x + g)), "no posterior draws"),
    "'g' produce several model-matrix columns"
  )
})

test_that("brms arguments are checked", {
  fit <- brms_fixture("uni")
  expect_error(relweights(fit, bootstrap = TRUE), "Bootstrap intervals")
  expect_error(relweights(fit, foo = 1), "Unknown argument")
  expect_error(relweights(fit, ndraws = -1), "ndraws")
  expect_error(relweights(fit, conf_level = 2), "conf_level")
  expect_error(relweights(fit, outcome = 1), "character vector")
})
