# Models and helpers shared by the tests.

# lavaan warns about some of these data sets (for example that the variances
# in mtcars differ by a factor of more than 1000); such warnings concern the
# fixtures, not relweights, so they are silenced when the models are fitted.
quiet_sem <- function(...) suppressWarnings(lavaan::sem(...))

hs_data <- lavaan::HolzingerSwineford1939

hs_model <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
  speed ~ visual + textual + ageyr + grade
"

pd_model <- "
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  dem60 ~ ind60
  dem65 ~ ind60 + dem60
  y1 ~~ y5
  y2 ~~ y4 + y6
  y3 ~~ y7
  y4 ~~ y8
  y6 ~~ y8
"

twolevel_model <- "
  level: 1
    fw =~ y1 + y2 + y3
    fw ~ x1 + x2 + x3
  level: 2
    fb =~ y1 + y2 + y3
    fb ~ w1 + w2
"

mtcars_predictors <- c("cyl", "disp", "hp", "gear")

# Data for a model in which the residual of y covaries with its predictor m
# (m is endogenous, with x as an instrument).
residual_cov_data <- function(n = 500) {
  set.seed(1)
  x <- stats::rnorm(n)
  u <- stats::rnorm(n)
  m <- x + u + stats::rnorm(n)
  y <- 0.5 * m + u + stats::rnorm(n)
  data.frame(x, m, y)
}
residual_cov_model <- "y ~ m\n m ~ x\n y ~~ m"

with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

layer_geoms <- function(p) {
  vapply(p$layers, function(l) class(l$geom)[1L], character(1))
}
