#' @rdname relweights
#' @export
relweights.lm <- function(object, bootstrap = FALSE, n_boot = 1000,
                          conf_level = 0.95, ...) {
  if (!identical(class(object), "lm")) {
    rw_stop("relweights() supports linear models of class 'lm' only; ",
            "objects of class ", rw_collapse(class(object)[1L]),
            " are not supported.")
  }
  rw_check_flag(bootstrap, "bootstrap")
  rw_check_boot_args(n_boot, conf_level)
  dots <- list(...)
  rw_check_dots(dots, c("parallel", "ncpus", "cl"), "linear models")

  dat <- rw_lm_data(object)
  if (length(dat$multi_column_terms)) {
    rw_warn("Term(s) ", rw_collapse(dat$multi_column_terms), " produce ",
            "several model-matrix columns, each of which receives its own ",
            "relative weight. For factors, polynomials and splines these ",
            "weights, and the weights of the other predictors, depend on the ",
            "contrasts or basis used to code the term.")
  }
  constant <- rw_constant_cols(dat$M)
  if (length(constant)) {
    rw_stop("Variable(s) ", rw_collapse(constant), " have zero variance in ",
            "the data used by the model.")
  }

  spec <- list(
    block = 1L, group = NA_character_, level = NA_character_,
    outcome = dat$outcome, predictors = dat$predictors, terms = dat$terms,
    n = nrow(dat$M), r2_model = summary(object)$r.squared,
    label = dat$outcome
  )
  fit <- rw_compute(rw_lm_cor(dat$M, dat$w), spec$outcome, spec$predictors,
                    label = spec$label)
  if (abs(fit$r2 - spec$r2_model) > 1e-8) {
    rw_warn("The R-squared implied by the correlations (",
            format(fit$r2, digits = 6), ") differs from the R-squared of the ",
            "model (", format(spec$r2_model, digits = 6), "); the ",
            "predictors may be close to collinear.")
  }

  info <- list(
    model_class = "lm", call = object$call, weighted = !is.null(dat$w),
    ngroups = 1L, nlevels = 1L, group_labels = character(),
    level_labels = character(), multi_column_terms = dat$multi_column_terms
  )
  out <- rw_build(list(spec), list(fit), info)
  if (bootstrap) {
    out <- rw_boot_lm(out, dat, spec, n_boot, conf_level, dots)
  }
  out
}

# Outcome, predictor columns and prior weights of a linear model.
rw_lm_data <- function(object) {
  tt <- stats::terms(object)
  if (attr(tt, "intercept") == 0L) {
    rw_stop("The model has no intercept. Relative weights are computed from ",
            "correlations, which assume a model with an intercept.")
  }
  mf <- stats::model.frame(object)
  if (!is.null(stats::model.offset(mf))) {
    rw_stop("The model has an offset, which relative weights cannot take ",
            "into account. Refit the model without the offset.")
  }
  cf <- stats::coef(object)
  if (anyNA(cf)) {
    rw_stop("The model has aliased coefficients: ",
            rw_collapse(names(cf)[is.na(cf)]), ". Remove the redundant ",
            "predictors and refit the model.")
  }
  X <- stats::model.matrix(object)
  assign <- attr(X, "assign")
  keep <- assign != 0L
  if (!any(keep)) {
    rw_stop("The model has no predictors.")
  }
  term_labels <- attr(tt, "term.labels")
  columns_per_term <- table(assign[keep])
  multi <- term_labels[as.integer(names(columns_per_term)[columns_per_term > 1L])]

  outcome <- names(mf)[attr(tt, "response")]
  predictors <- colnames(X)[keep]
  y <- as.numeric(stats::model.response(mf, "numeric"))
  M <- cbind(y, X[, keep, drop = FALSE])
  dimnames(M) <- list(NULL, c(outcome, predictors))

  w <- stats::model.weights(mf)
  if (!is.null(w)) {
    positive <- w > 0
    M <- M[positive, , drop = FALSE]
    w <- as.numeric(w[positive])
  }
  list(outcome = outcome, predictors = predictors,
       terms = term_labels[assign[keep]], multi_column_terms = multi,
       M = M, w = w)
}

rw_constant_cols <- function(M) {
  constant <- vapply(seq_len(ncol(M)), function(j) {
    v <- M[, j]
    all(v == v[1L])
  }, logical(1))
  colnames(M)[constant]
}

rw_lm_cor <- function(M, w = NULL) {
  if (is.null(w)) {
    stats::cor(M)
  } else {
    stats::cov.wt(M, wt = w, cor = TRUE, method = "ML")$cor
  }
}

rw_boot_lm <- function(x, dat, spec, n_boot, conf_level, dots) {
  n_stat <- 2L * length(spec$predictors) + 1L
  k <- ncol(dat$M)
  weighted <- !is.null(dat$w)
  data <- if (weighted) cbind(dat$M, dat$w) else dat$M

  statistic <- function(data, indices) {
    d <- data[indices, , drop = FALSE]
    M <- d[, seq_len(k), drop = FALSE]
    w <- if (weighted) d[, k + 1L] else NULL
    if (length(rw_constant_cols(M))) {
      return(rep(NA_real_, n_stat))
    }
    fit <- tryCatch(
      suppressWarnings(rw_compute(rw_lm_cor(M, w), spec$outcome,
                                  spec$predictors)),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(rep(NA_real_, n_stat))
    }
    c(fit$raw, fit$rescaled, fit$r2)
  }

  args <- c(list(data = data, statistic = statistic, R = n_boot), dots)
  b <- do.call(boot::boot, args)
  rw_add_intervals(x, b$t, b$t0, rw_stat_index(list(spec)), conf_level,
                   method = "boot::boot")
}
