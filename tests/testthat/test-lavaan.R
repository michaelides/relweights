test_that("latent regressions reproduce lavaan's R2 and standardized coefficients", {
  fit <- quiet_sem(hs_model, data = hs_data)
  rw <- relweights(fit)
  std <- lavaan::standardizedSolution(fit)
  std <- std[std$op == "~", ]

  expect_s3_class(rw, "relweights")
  expect_equal(rw$weights$predictor, c("visual", "textual", "ageyr", "grade"))
  expect_equal(rw$weights$term, rw$weights$predictor)
  expect_equal(rw$weights$beta, std$est.std)
  expect_equal(rw$r2$r2, unname(lavaan::lavInspect(fit, "rsquare")["speed"]))
  expect_equal(rw$r2$r2_model, rw$r2$r2)
  expect_equal(rw$r2$n, lavaan::lavInspect(fit, "nobs"))
  expect_true(all(is.na(rw$weights$group)))
  expect_true(all(is.na(rw$weights$level)))

  R <- unclass(lavaan::lavInspect(fit, "cor.all"))
  expect_equal(rw$weights$r, unname(R[rw$weights$predictor, "speed"]))
  expect_identical(rw$info$model_class, "lavaan")
})

test_that("every regression is analysed and `outcome` selects some", {
  fit <- quiet_sem(pd_model, data = lavaan::PoliticalDemocracy)
  rw <- relweights(fit)
  rsquare <- lavaan::lavInspect(fit, "rsquare")

  expect_equal(rw$r2$outcome, c("dem60", "dem65"))
  expect_equal(rw$r2$r2, unname(rsquare[c("dem60", "dem65")]))
  expect_equal(rw$weights$outcome, c("dem60", "dem65", "dem65"))
  expect_equal(rw$weights$rescaled[1], 100)
  expect_named(rw$matrices, c("dem60", "dem65"))

  only <- relweights(fit, outcome = "dem65")
  expect_equal(only$r2$outcome, "dem65")
  expect_equal(only$weights$raw, rw$weights$raw[2:3])

  expect_error(relweights(fit, outcome = "dem70"),
               "'dem70' are not regressed.*'dem60', 'dem65'")
  expect_error(relweights(fit, outcome = 1), "character vector")
})

test_that("groups are analysed separately and labelled", {
  fit <- quiet_sem(hs_model, data = hs_data, group = "school")
  rw <- relweights(fit)
  rsquare <- lavaan::lavInspect(fit, "rsquare")

  expect_equal(rw$r2$group, c("Pasteur", "Grant-White"))
  expect_equal(rw$r2$block, 1:2)
  expect_equal(rw$r2$n, lavaan::lavInspect(fit, "nobs"))
  expect_equal(rw$r2$r2, unname(c(rsquare$Pasteur["speed"],
                                  rsquare$`Grant-White`["speed"])))
  expect_equal(rw$weights$group, rep(c("Pasteur", "Grant-White"), each = 4))
  expect_named(rw$matrices, c("speed [Pasteur]", "speed [Grant-White]"))
  expect_equal(rw$info$ngroups, 2L)
  expect_equal(rw$info$group_labels, c("Pasteur", "Grant-White"))
})

test_that("levels are analysed separately with the right sample sizes", {
  fit <- quiet_sem(twolevel_model, data = lavaan::Demo.twolevel,
                   cluster = "cluster")
  rw <- relweights(fit)
  rsquare <- lavaan::lavInspect(fit, "rsquare")

  expect_equal(rw$r2$level, c("within", "cluster"))
  expect_equal(rw$r2$outcome, c("fw", "fb"))
  expect_equal(rw$r2$n, c(2500L, 200L))
  expect_equal(rw$r2$r2, unname(c(rsquare$within["fw"], rsquare$cluster["fb"])))
  expect_equal(rw$weights$predictor, c("x1", "x2", "x3", "w1", "w2"))
  expect_equal(rw$info$nlevels, 2L)
})

test_that("blocks of multigroup two-level models get group and level labels", {
  d <- lavaan::Demo.twolevel
  d$g <- ifelse(d$cluster %% 2 == 0, "even", "odd")
  block <- "
    level: 1
      y1 ~ x1 + x2
    level: 2
      y1 ~ w1 + w2
  "
  model <- paste("group: odd", block, "group: even", block, sep = "\n")
  fit <- quiet_sem(model, data = d, cluster = "cluster", group = "g")
  rw <- relweights(fit)
  rsquare <- unlist(lapply(lavaan::lavInspect(fit, "rsquare"),
                           function(r) r["y1"]))

  expect_equal(rw$r2$group, rep(c("odd", "even"), each = 2))
  expect_equal(rw$r2$level, rep(c("within", "cluster"), 2))
  expect_equal(rw$r2$n, c(1000L, 100L, 1500L, 100L))
  expect_equal(rw$r2$r2, unname(rsquare))
  expect_equal(rw$info$group_labels, c("odd", "even"))
  expect_equal(rw$info$level_labels, c("within", "cluster"))
})

test_that("models fitted to summary statistics are supported", {
  S <- cov(mtcars[, c("mpg", mtcars_predictors)])
  fit <- quiet_sem("mpg ~ cyl + disp + hp + gear", sample.cov = S,
                   sample.nobs = 32)
  rw <- relweights(fit)
  reference <- relweights(lm(mpg ~ cyl + disp + hp + gear, data = mtcars))
  expect_equal(rw$weights$raw, reference$weights$raw, tolerance = 1e-6)
  expect_equal(rw$r2$n, 32L)
})

test_that("ordered indicators are analysed on the latent-response scale", {
  d <- hs_data
  for (v in paste0("x", 1:9)) d[[v]] <- cut(d[[v]], 3, labels = FALSE)
  model <- "
    visual  =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
    speed ~ visual + textual + ageyr
  "
  conditional <- quiet_sem(model, data = d, ordered = paste0("x", 1:9))
  expect_error(relweights(conditional), "conditional.x = FALSE")

  fit <- quiet_sem(model, data = d, ordered = paste0("x", 1:9),
                   conditional.x = FALSE)
  rw <- relweights(fit)
  expect_equal(rw$r2$r2, unname(lavaan::lavInspect(fit, "rsquare")["speed"]))
})

test_that("a residual covariance between outcome and predictor is flagged", {
  fit <- quiet_sem(residual_cov_model, data = residual_cov_data())

  expect_warning(rw <- relweights(fit), "For 'y'.*residual")
  r2 <- rw$r2[rw$r2$outcome == "y", ]
  expect_gt(abs(r2$r2 - r2$r2_model), 0.05)
  expect_equal(r2$r2_model, unname(lavaan::lavInspect(fit, "rsquare")["y"]))
  R <- unclass(lavaan::lavInspect(fit, "cor.all"))
  expect_equal(r2$r2, R["y", "m"]^2)
})

test_that("problematic lavaan models are flagged", {
  not_converged <- quiet_sem(
    "visual =~ x1 + x2 + x3\n textual =~ x4 + x5 + x6\n textual ~ visual",
    data = hs_data, control = list(iter.max = 2)
  )
  expect_error(relweights(not_converged), "has not converged")

  inadmissible <- quiet_sem(
    "textual =~ x4 + x5 + x6\n speed =~ x7 + x8 + x9\n x5 ~~ -0.05*x5\n speed ~ textual + ageyr",
    data = hs_data
  )
  expect_warning(rw <- relweights(inadmissible), "not admissible")
  expect_equal(rw$r2$r2,
               unname(lavaan::lavInspect(inadmissible, "rsquare")["speed"]))

  measurement_only <- quiet_sem("visual =~ x1 + x2 + x3", data = hs_data)
  expect_error(relweights(measurement_only), "no regressions")
})

test_that("unknown arguments are errors for lavaan models", {
  fit <- quiet_sem("x9 ~ x7 + x8", data = hs_data)
  expect_error(relweights(fit, foo = 1), "lavaan models: 'foo'")
})
