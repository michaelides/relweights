named_cor <- function(values, names) {
  dimnames(values) <- list(names, names)
  values
}

test_that("orthogonal predictors receive their squared correlations", {
  R <- named_cor(diag(3), c("y", "a", "b"))
  R["y", "a"] <- R["a", "y"] <- 0.5
  R["y", "b"] <- R["b", "y"] <- 0.3

  fit <- rw_compute(R, "y", c("a", "b"))

  expect_equal(unname(fit$raw), c(0.25, 0.09))
  expect_equal(fit$r2, 0.34)
  expect_equal(unname(fit$beta), c(0.5, 0.3))
  expect_equal(unname(fit$orth), c(0.5, 0.3))
  expect_equal(unname(fit$rescaled), 100 * c(0.25, 0.09) / 0.34)
})

test_that("raw weights sum to R2 and rescaled weights to 100", {
  R <- cor(mtcars[, c("mpg", mtcars_predictors)])
  fit <- rw_compute(R, "mpg", mtcars_predictors)
  beta <- solve(R[-1, -1], R[-1, 1])

  expect_equal(sum(fit$raw), fit$r2)
  expect_equal(sum(fit$rescaled), 100)
  expect_true(all(fit$raw >= 0))
  expect_equal(fit$beta, beta)
  expect_equal(fit$r2, sum(beta * R[-1, 1]))
  expect_equal(fit$r, R[-1, 1])
  # Lambda reproduces the predictor correlations: Lambda %*% Lambda = Rxx.
  expect_equal(fit$lambda %*% fit$lambda, R[-1, -1])
  expect_named(fit$raw, mtcars_predictors)
})

test_that("the order of the variables in the matrix does not matter", {
  R <- cor(mtcars[, c(mtcars_predictors, "mpg")])
  fit <- rw_compute(R, "mpg", mtcars_predictors)
  reference <- rw_compute(cor(mtcars[, c("mpg", mtcars_predictors)]), "mpg",
                          mtcars_predictors)
  expect_equal(fit$raw, reference$raw)
  expect_equal(rownames(fit$R), c("mpg", mtcars_predictors))
})

test_that("weights do not change when variables are rescaled", {
  rescaled <- transform(mtcars, disp = disp / 100 + 3, mpg = -2 * mpg)
  original <- relweights(lm(mpg ~ cyl + disp + hp, data = mtcars))
  changed <- relweights(lm(mpg ~ cyl + disp + hp, data = rescaled))

  expect_equal(changed$weights$raw, original$weights$raw)
  expect_equal(changed$weights$sign_beta, -original$weights$sign_beta)
  expect_equal(changed$weights$sign_orth, -original$weights$sign_orth)
})

test_that("collinear predictors are an error", {
  R <- named_cor(matrix(c(1, 0.5, 0.5,
                          0.5, 1, 1,
                          0.5, 1, 1), 3), c("y", "a", "b"))
  expect_error(rw_compute(R, "y", c("a", "b")), "collinear")
})

test_that("invalid correlation matrices are errors", {
  not_psd <- named_cor(matrix(c(1, 0.9, -0.9,
                                0.9, 1, 0.9,
                                -0.9, 0.9, 1), 3), c("y", "a", "b"))
  expect_error(rw_compute(not_psd, "y", c("a", "b")),
               "not positive semi-definite")

  with_na <- named_cor(diag(2), c("y", "a"))
  with_na["y", "a"] <- with_na["a", "y"] <- NA
  expect_error(rw_compute(with_na, "y", "a"), "non-finite")

  expect_error(rw_compute(named_cor(diag(2), c("y", "a")), "y", "b"),
               "missing from the correlation matrix")
})

test_that("a zero R2 leaves rescaled weights undefined", {
  R <- named_cor(diag(3), c("y", "a", "b"))
  expect_warning(fit <- rw_compute(R, "y", c("a", "b")), "R-squared")
  expect_equal(fit$r2, 0)
  expect_true(all(is.na(fit$rescaled)))
})

test_that("percentile intervals skip failed replicates", {
  expect_equal(rw_percentile(c(1:99, NA, NaN), 0.9),
               unname(stats::quantile(1:99, c(0.05, 0.95), type = 6)))
  expect_equal(rw_percentile(c(NA_real_, NA_real_), 0.95),
               c(NA_real_, NA_real_))
})

test_that("labels combine outcome, group and level", {
  expect_equal(rw_label("y"), "y")
  expect_equal(rw_label(c("y", "z")), c("y", "z"))
  expect_equal(rw_label(c("y", "y"), c("g1", "g2"), NA_character_),
               c("y [g1]", "y [g2]"))
  expect_equal(rw_label("y", "g1", "within"), "y [g1, within]")
  expect_equal(rw_label("y", NA_character_, "cluster"), "y [cluster]")
})
