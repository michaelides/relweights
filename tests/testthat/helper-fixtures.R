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

# Bayesian models ----------------------------------------------------------

# Small brms models saved by data-raw/brms-fixtures.R (fitting them would
# compile Stan programs).
brms_fixture <- function(name) {
  readRDS(test_path("fixtures", paste0("brms-", name, ".rds")))
}

# brms models without draws, for testing which models are accepted; they are
# built without compiling anything.
brms_empty_data <- data.frame(
  x = sin(1:40), z = cos(1.3 * (1:40)), g = rep(letters[1:4], 10),
  w = rep(c(0.5, 1, 1.5, 2), 10), o = rep(1:3, length.out = 40),
  b = rep(0:1, 20)
)
brms_empty_data$y <- brms_empty_data$x + 0.5 * brms_empty_data$z +
  sin(2.1 * (1:40))

empty_brm <- function(formula, family = gaussian()) {
  suppressMessages(suppressWarnings(
    brms::brm(formula, data = brms_empty_data, family = family, empty = TRUE)
  ))
}

# blavaan models use precompiled Stan programs and take a few seconds.
quiet_bsem <- function(..., burnin = 300, sample = 300) {
  # blavaan::bsem() calls blavaan() in the frame it is called from, so the
  # function must be visible there even when blavaan is not attached.
  blavaan <- blavaan::blavaan
  out <- NULL
  utils::capture.output(out <- suppressMessages(suppressWarnings(
    blavaan::bsem(..., burnin = burnin, sample = sample, n.chains = 2,
                  seed = 1, bcontrol = list(cores = 1, refresh = 0))
  )))
  out
}

with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

layer_geoms <- function(p) {
  vapply(p$layers, function(l) class(l$geom)[1L], character(1))
}
