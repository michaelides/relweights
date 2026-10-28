skip_if_not_installed("blavaan")
skip_if_not_installed("posterior")

latent_model <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
  speed ~ visual + textual + ageyr
"

test_that("at the point estimates the per-draw code reproduces lavaan", {
  fit <- quiet_bsem(latent_model, data = hs_data)
  specs <- rw_lavaan_specs(fit)
  res <- rw_blavaan_draw(fit, specs, lavaan::coef(fit),
                         rw_blavaan_resvar(fit, specs))
  ref <- rw_fit_specs(specs, rw_lavaan_cors(fit))

  expect_equal(res$cors[[1]], unclass(lavaan::lavInspect(fit, "cor.all")),
               tolerance = 1e-10)
  expect_equal(res$fits[[1]]$raw, ref[[1]]$raw, tolerance = 1e-10)
  expect_equal(res$fits[[1]]$beta, ref[[1]]$beta, tolerance = 1e-10)
  expect_equal(res$r2_model,
               unname(lavaan::lavInspect(fit, "rsquare")["speed"]),
               tolerance = 1e-10)
})

test_that("latent regressions are analysed in every draw", {
  fit <- quiet_bsem(latent_model, data = hs_data)
  # Short chains of latent-variable models occasionally get stuck (for
  # example where a factor variance is close to zero); relweights() then
  # warns, and the comparison with maximum likelihood is skipped below.
  converged <- max(blavaan::blavInspect(fit, "psrf")) < 1.05
  rw <- if (converged) relweights(fit) else suppressWarnings(relweights(fit))
  draws <- rw$draws$regressions$speed

  expect_identical(rw$info$model_class, "blavaan")
  expect_equal(rw$weights$predictor, c("visual", "textual", "ageyr"))
  expect_equal(rw$r2$n, 301L)
  expect_equal(nrow(draws$raw), 600L)
  expect_equal(draws$r2, draws$r2_model, tolerance = 1e-8)
  expect_equal(rowSums(draws$raw), draws$r2)
  expect_equal(sum(rw$weights$rescaled), 100)
  expect_true(all(is.finite(rw$weights$rhat)))
  expect_true(all(rw$weights$ess_bulk > 0))

  thinned <- suppressWarnings(relweights(fit, ndraws = 100))
  expect_equal(thinned$info$ndraws, 100L)
  expect_equal(unique(thinned$draws$chain), 1:2)
  expect_equal(thinned$draws$regressions$speed$raw,
               draws$raw[c(rw_thin_iterations(300, 2, 100),
                           300 + rw_thin_iterations(300, 2, 100)), ])

  out <- capture.output(print(rw))
  expect_true(any(grepl("Model: blavaan model", out, fixed = TRUE)))
  expect_true(any(grepl("rescaled CrI", out, fixed = TRUE)))
  expect_error(relweights(fit, bootstrap = TRUE), "Bootstrap intervals")

  skip_if_not(converged, "the MCMC chains did not converge")
  ref <- relweights(lavaan::sem(latent_model, data = hs_data))
  expect_lt(max(abs(rw$weights$rescaled - ref$weights$rescaled)), 5)
  expect_lt(abs(rw$r2$r2 - ref$r2$r2), 0.05)
  expect_true(all(rw$weights$rhat < 1.05))
})

test_that("groups and levels are analysed separately", {
  fit <- quiet_bsem("x9 ~ x7 + x8 + ageyr", data = hs_data, group = "school",
                    burnin = 200, sample = 200)
  rw <- relweights(fit)
  expect_equal(rw$r2$group, c("Pasteur", "Grant-White"))
  expect_equal(rw$r2$n, c(156L, 145L))
  expect_equal(rw$weights$group, rep(c("Pasteur", "Grant-White"), each = 3))
  expect_named(rw$draws$regressions, c("x9 [Pasteur]", "x9 [Grant-White]"))
  expect_equal(rw$draws$regressions[[2]]$r2,
               rw$draws$regressions[[2]]$r2_model, tolerance = 1e-8)
  p <- autoplot(rw, type = "groups")
  expect_no_warning(ggplot2::ggplot_build(p))
  expect_true("GeomPointrange" %in% layer_geoms(p))

  d <- lavaan::Demo.twolevel
  d <- d[d$cluster <= 60, ]
  model <- "
    level: 1
      y1 ~ x1 + x2
    level: 2
      y1 ~ w1 + w2
  "
  fit <- quiet_bsem(model, data = d, cluster = "cluster")
  rw <- relweights(fit)
  expect_equal(rw$r2$level, c("within", "cluster"))
  expect_equal(rw$r2$n, c(750L, 60L))
  expect_equal(rw$info$nlevels, 2L)

  specs <- rw_lavaan_specs(fit)
  res <- rw_blavaan_draw(fit, specs, lavaan::coef(fit),
                         rw_blavaan_resvar(fit, specs))
  expect_equal(res$cors, lapply(lavaan::lavInspect(fit, "cor.all"), unclass),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(res$r2_model,
               unname(vapply(lavaan::lavInspect(fit, "rsquare"),
                             function(r) r[["y1"]], numeric(1))),
               tolerance = 1e-10)
})

test_that("an outcome whose residual covaries with a predictor is flagged", {
  data <- residual_cov_data(300)
  fit <- quiet_bsem(residual_cov_model, data = data)
  expect_warning(rw <- relweights(fit), "For 'y'.*residual")
  y <- rw$r2[rw$r2$outcome == "y", ]
  expect_gt(y$r2 - y$r2_model, 0.05)
  expect_equal(relweights(fit, outcome = "m")$r2$outcome, "m")
})
