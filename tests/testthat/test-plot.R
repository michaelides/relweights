lm_fit <- lm(mpg ~ cyl + disp + hp + wt + qsec, data = mtcars)

test_that("every plot type, scale and sign can be built", {
  rw <- relweights(lm_fit)
  for (type in c("weights", "r2", "compare")) {
    for (scale in c("rescaled", "raw")) {
      for (sign in c("none", "beta", "orthogonal")) {
        p <- autoplot(rw, type = type, scale = scale, sign = sign)
        expect_s3_class(p, "ggplot")
        expect_no_warning(ggplot2::ggplot_build(p))
      }
    }
  }
})

test_that("signed weights take the sign of the chosen coefficient", {
  rw <- relweights(lm(mpg ~ ., data = mtcars))
  w <- rw$weights

  by_beta <- rw_weight_values(rw, "rescaled", "beta")
  expect_equal(by_beta$value, w$rescaled * w$sign_beta)
  by_orth <- rw_weight_values(rw, "raw", "orthogonal")
  expect_equal(by_orth$value, w$raw * w$sign_orth)
  expect_equal(rw_weight_values(rw, "raw", "none")$value, w$raw)

  p <- autoplot(rw, sign = "beta")
  bars <- ggplot2::layer_data(p, which(layer_geoms(p) == "GeomCol"))
  # disp has a positive standardized coefficient but a negative
  # orthogonal-counterpart coefficient.
  disp <- by_beta$value[w$predictor == "disp"]
  expect_gt(disp, 0)
  expect_equal(sum(bars$xmin < 0), sum(w$sign_beta < 0))
  expect_equal(sum(rw_weight_values(rw, "raw", "orthogonal")$value < 0),
               sum(w$sign_orth < 0))
})

test_that("signed intervals are mirrored for negative weights", {
  set.seed(5)
  rw <- relweights(lm(mpg ~ wt + qsec, data = mtcars), bootstrap = TRUE,
                   n_boot = 50)
  d <- rw_weight_values(rw, "raw", "beta")
  negative <- d$sign_beta < 0
  expect_true(any(negative))
  expect_equal(d$lower[negative], -rw$weights$raw_upper[negative])
  expect_equal(d$upper[negative], -rw$weights$raw_lower[negative])
  expect_true(all(d$lower <= d$value & d$value <= d$upper))
})

test_that("error bars are drawn only for bootstrap results", {
  rw <- relweights(lm_fit)
  expect_false("GeomErrorbar" %in% layer_geoms(autoplot(rw)))
  expect_false("GeomErrorbar" %in% layer_geoms(autoplot(rw, type = "r2")))

  set.seed(6)
  rwb <- relweights(lm_fit, bootstrap = TRUE, n_boot = 50)
  expect_true("GeomErrorbar" %in% layer_geoms(autoplot(rwb)))
  expect_true("GeomErrorbar" %in% layer_geoms(autoplot(rwb, type = "r2")))
  expect_false("GeomErrorbar" %in% layer_geoms(autoplot(rwb, ci = FALSE)))
  expect_no_warning(ggplot2::ggplot_build(autoplot(rwb, sign = "orthogonal")))
  expect_no_warning(ggplot2::ggplot_build(autoplot(rwb, type = "r2")))
})

test_that("bars are sorted by weight", {
  rw <- relweights(lm_fit)
  p <- autoplot(rw)
  order_on_axis <- rw_strip_within(levels(p$data$y))
  expected <- rw$weights$predictor[order(rw$weights$rescaled)]
  expect_equal(order_on_axis, expected)
})

test_that("the R2 plot stacks the raw weights up to R2", {
  rw <- relweights(lm_fit)
  p <- autoplot(rw, type = "r2")
  tiles <- ggplot2::layer_data(p, which(layer_geoms(p) == "GeomTile"))
  expect_equal(max(tiles$xmax), rw$r2$r2)
  expect_equal(tiles$xmax - tiles$xmin, rw$weights$raw)
})

test_that("the comparison plot shows three indices that each sum to 100", {
  rw <- relweights(lm_fit)
  p <- autoplot(rw, type = "compare")
  totals <- tapply(p$data$value, p$data$measure, sum)
  expect_equal(as.vector(totals), rep(100, 3))
  expect_setequal(levels(p$data$measure),
                  c("Relative weight", "Squared standardized coefficient",
                    "Squared correlation"))
})

test_that("lavaan results are plotted by regression, group and level", {
  fit <- quiet_sem(hs_model, data = hs_data, group = "school")
  rw <- relweights(fit)
  for (type in c("weights", "r2", "compare", "groups")) {
    p <- autoplot(rw, type = type, sign = "beta")
    expect_s3_class(p, "ggplot")
    expect_no_warning(ggplot2::ggplot_build(p))
  }
  p <- autoplot(rw, type = "groups")
  expect_equal(levels(p$data$block_label), c("Grant-White", "Pasteur"))
  expect_equal(nlevels(autoplot(rw)$data$panel), 2L)

  multilevel <- relweights(quiet_sem(twolevel_model,
                                     data = lavaan::Demo.twolevel,
                                     cluster = "cluster"))
  p <- autoplot(multilevel, type = "groups")
  expect_identical(p$scales$get_scales("colour")$name, "Level")
  expect_no_warning(ggplot2::ggplot_build(p))
})

test_that("plots can be restricted to some outcomes", {
  fit <- quiet_sem(pd_model, data = lavaan::PoliticalDemocracy)
  rw <- relweights(fit)
  p <- autoplot(rw, outcome = "dem65")
  expect_equal(unique(p$data$outcome), "dem65")
  expect_error(autoplot(rw, outcome = "dem70"), "'dem70' are not part")
})

test_that("invalid plot requests are informative errors", {
  rw <- relweights(lm_fit)
  expect_error(autoplot(rw, type = "groups"), "single group and level")
  expect_error(autoplot(rw, type = "pie"), "should be one of")
  expect_error(autoplot(rw, colour = "red"), "Unused argument\\(s\\).*colour")
  expect_error(autoplot(rw, ci = NA), "TRUE or FALSE")
})

test_that("plot() draws the plot and returns it invisibly", {
  rw <- relweights(lm_fit)
  with_null_device({
    expect_invisible(p <- plot(rw, type = "compare"))
  })
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Importance indices")
})
