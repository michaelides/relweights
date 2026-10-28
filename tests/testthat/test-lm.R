test_that("relweights() reproduces the R2 and the standardized coefficients of lm()", {
  fit <- lm(mpg ~ cyl + disp + hp + gear, data = mtcars)
  rw <- relweights(fit)
  scaled <- lm(scale(mpg) ~ scale(cyl) + scale(disp) + scale(hp) + scale(gear),
               data = mtcars)

  expect_s3_class(rw, "relweights")
  expect_equal(rw$r2$r2, summary(fit)$r.squared)
  expect_equal(rw$r2$r2_model, summary(fit)$r.squared)
  expect_equal(rw$r2$n, 32L)
  expect_equal(rw$weights$beta, unname(coef(scaled)[-1]))
  expect_equal(rw$weights$r, unname(cor(mtcars[mtcars_predictors], mtcars$mpg)[, 1]))
  expect_equal(rw$weights$predictor, names(coef(fit))[-1])
  expect_equal(rw$weights$sign_beta, as.integer(sign(coef(fit)[-1])),
               ignore_attr = TRUE)
  expect_equal(sum(rw$weights$raw), rw$r2$r2)
  expect_null(rw$boot)
  expect_identical(rw$info$model_class, "lm")
  expect_false(rw$info$weighted)
  expect_named(rw$matrices, "mpg")
  expect_equal(rownames(rw$matrices$mpg$R), c("mpg", mtcars_predictors))
})

test_that("prior weights give weighted correlations and drop zero weights", {
  w <- rep(c(0, 1, 2, 3), length.out = nrow(mtcars))
  fit <- lm(mpg ~ wt + hp + qsec, data = mtcars, weights = w)
  rw <- relweights(fit)

  expect_equal(rw$r2$r2, summary(fit)$r.squared)
  expect_equal(rw$r2$n, sum(w > 0))
  expect_true(rw$info$weighted)

  positive <- w > 0
  refit <- lm(mpg ~ wt + hp + qsec, data = mtcars[positive, ],
              weights = w[positive])
  expect_equal(relweights(refit)$weights$raw, rw$weights$raw)
})

test_that("rows dropped by the model are not used", {
  d <- mtcars
  d$wt[1:3] <- NA
  omit <- relweights(lm(mpg ~ wt + hp, data = d))
  exclude <- relweights(lm(mpg ~ wt + hp, data = d, na.action = na.exclude))
  complete <- relweights(lm(mpg ~ wt + hp, data = d[-(1:3), ]))

  expect_equal(omit$r2$n, 29L)
  expect_equal(omit$weights$raw, complete$weights$raw)
  expect_equal(exclude$weights$raw, complete$weights$raw)
})

test_that("transformed responses and predictors keep their model names", {
  rw <- relweights(lm(log(mpg) ~ log(hp) + I(wt^2), data = mtcars))
  expect_equal(rw$r2$outcome, "log(mpg)")
  expect_equal(rw$weights$predictor, c("log(hp)", "I(wt^2)"))
  expect_equal(rw$weights$term, c("log(hp)", "I(wt^2)"))
})

test_that("terms with several columns get one weight per column and a warning", {
  expect_warning(
    rw <- relweights(lm(mpg ~ factor(cyl) + wt + hp, data = mtcars)),
    "factor\\(cyl\\).*contrasts or basis"
  )
  expect_equal(rw$weights$predictor, c("factor(cyl)6", "factor(cyl)8", "wt", "hp"))
  expect_equal(rw$weights$term, c("factor(cyl)", "factor(cyl)", "wt", "hp"))
  expect_equal(rw$info$multi_column_terms, "factor(cyl)")
  expect_equal(rw$r2$r2, summary(lm(mpg ~ factor(cyl) + wt + hp, data = mtcars))$r.squared)

  expect_warning(relweights(lm(mpg ~ poly(wt, 2) + hp, data = mtcars)),
                 "poly\\(wt, 2\\)")
})

test_that("two-level factors and interactions give one column without warning", {
  expect_no_warning(rw <- relweights(lm(mpg ~ factor(am) + wt, data = mtcars)))
  expect_equal(rw$weights$predictor, c("factor(am)1", "wt"))
  expect_length(rw$info$multi_column_terms, 0L)

  expect_no_warning(rw <- relweights(lm(mpg ~ wt * hp, data = mtcars)))
  expect_equal(rw$weights$term, c("wt", "hp", "wt:hp"))
})

test_that("the reference level of a two-level factor only flips signs", {
  d <- mtcars
  d$am_f <- factor(d$am)
  first <- relweights(lm(mpg ~ am_f + wt, data = d))
  d$am_f <- relevel(d$am_f, ref = "1")
  second <- relweights(lm(mpg ~ am_f + wt, data = d))
  expect_equal(second$weights$raw, first$weights$raw)
  expect_equal(second$weights$sign_beta[1], -first$weights$sign_beta[1])
})

test_that("unsupported linear models are informative errors", {
  expect_error(relweights(lm(mpg ~ 0 + wt + hp, data = mtcars)), "no intercept")
  expect_error(relweights(lm(mpg ~ wt + offset(hp), data = mtcars)), "offset")
  expect_error(relweights(lm(mpg ~ wt, offset = hp, data = mtcars)), "offset")

  d <- mtcars
  d$wt2 <- 2 * d$wt
  expect_error(relweights(lm(mpg ~ wt + wt2 + hp, data = d)),
               "aliased coefficients: 'wt2'")
  expect_error(relweights(lm(mpg ~ 1, data = mtcars)), "no predictors")

  fit <- lm(mpg ~ wt + hp, data = mtcars)
  class(fit) <- c("other", "lm")
  expect_error(relweights(fit), "class 'other' are not supported")
})

test_that("objects without a method are informative errors", {
  expect_error(relweights(mtcars), "class 'data.frame' are not supported")
  expect_error(relweights(1:3), "class 'integer' are not supported")
})

test_that("unknown arguments are errors", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_error(relweights(fit, foo = 1), "Unknown argument\\(s\\).*'foo'")
  expect_error(relweights(fit, outcome = "mpg"), "'outcome'")
  expect_error(relweights(fit, FALSE, 1000, 0.95, "multicore"),
               "must be named")
  expect_error(relweights(fit, bootstrap = 1), "TRUE or FALSE")
})
