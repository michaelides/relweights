test_that("print() shows R2, the weights and the signs", {
  rw <- relweights(lm(mpg ~ cyl + disp + hp + gear, data = mtcars))
  out <- capture.output(res <- print(rw))

  expect_identical(res, rw)
  expect_match(out[1], "Relative weight analysis \\(Johnson, 2000\\)")
  expect_true(any(grepl("Model: lm(formula = mpg ~ cyl + disp + hp + gear",
                        out, fixed = TRUE)))
  expect_true(any(grepl("R-squared = 0.779, n = 32", out, fixed = TRUE)))
  expect_true(any(grepl("^\\s+cyl 0\\.228\\s+29\\.3 -0\\.257 -0\\.852\\s+-/-$", out)))
  expect_false(any(grepl("CI", out)))
  expect_false(any(grepl("differ", out)))
})

test_that("print() flags predictors whose two signs differ", {
  rw <- relweights(lm(mpg ~ ., data = mtcars))
  expect_output(print(rw), "The two signs differ for: disp \\(mpg\\)")
  expect_output(print(rw), "\\+/-")
})

test_that("print() shows bootstrap intervals", {
  set.seed(1)
  rw <- relweights(lm(mpg ~ wt + hp, data = mtcars), bootstrap = TRUE,
                   n_boot = 50)
  out <- capture.output(print(rw))
  expect_true(any(grepl("raw CI", out)))
  expect_true(any(grepl("rescaled CI", out)))
  expect_true(any(grepl("R-squared = 0.827 [", out, fixed = TRUE)))
  expect_true(any(grepl("95% percentile bootstrap intervals from 50 replicates",
                        out, fixed = TRUE)))
})

test_that("print() describes lavaan models by group and level", {
  rw <- relweights(quiet_sem(hs_model, data = hs_data, group = "school"))
  out <- capture.output(print(rw))
  expect_true(any(grepl("Model: lavaan model; groups: Pasteur, Grant-White",
                        out, fixed = TRUE)))
  expect_true(any(grepl("Outcome: speed [Grant-White]", out, fixed = TRUE)))

  fit <- quiet_sem(residual_cov_model, data = residual_cov_data())
  rw <- suppressWarnings(relweights(fit))
  expect_output(print(rw), "R-squared = 0\\.534 \\(model: 0\\.443\\)")
})

test_that("as.data.frame() returns the weights or the R2 table", {
  rw <- relweights(lm(mpg ~ wt + hp, data = mtcars))
  expect_identical(as.data.frame(rw), rw$weights)
  expect_identical(as.data.frame(rw, what = "r2"), rw$r2)
  expect_equal(rownames(as.data.frame(rw, row.names = c("a", "b"))),
               c("a", "b"))
  expect_error(as.data.frame(rw, what = "matrices"), "should be one of")
  expect_named(rw$weights, c("block", "group", "level", "outcome", "predictor",
                             "term", "raw", "rescaled", "beta", "orth", "r",
                             "sign_beta", "sign_orth"))
  expect_named(rw$r2, c("block", "group", "level", "outcome", "n", "r2",
                        "r2_model"))
})
