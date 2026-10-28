# Internal computational core shared by the lm and lavaan methods.
#
# A "regression specification" (spec) is a list with elements
#   block, group, level, outcome, predictors, terms, n, r2_model, label
# and a "fit" is the result of rw_compute() for one spec.

rw_stop <- function(...) stop(paste0(...), call. = FALSE)

rw_warn <- function(...) warning(paste0(...), call. = FALSE)

# Johnson's (2000) relative weights from a correlation matrix.
#
# With Rxx = V diag(d) V' (eigendecomposition of the predictor correlations),
# Lambda = V diag(sqrt(d)) V' holds the correlations between the predictors
# and their orthogonal counterparts Z. The outcome is regressed on Z
# (orth = Lambda^-1 r_xy), and the variance each Z_k explains is shared out
# among the predictors in proportion to Lambda_jk^2.
rw_compute <- function(R, outcome, predictors, label = outcome, tol = 1e-10) {
  vars <- c(outcome, predictors)
  R <- unclass(R)
  missing_vars <- setdiff(vars, rownames(R))
  if (length(missing_vars)) {
    rw_stop("Variable(s) ", rw_collapse(missing_vars),
            " are missing from the correlation matrix of '", label, "'.")
  }
  R <- R[vars, vars, drop = FALSE]
  storage.mode(R) <- "double"
  if (any(!is.finite(R))) {
    rw_stop("The correlation matrix of '", label,
            "' and its predictors contains missing or non-finite values.")
  }
  R <- (R + t(R)) / 2
  Rxx <- R[predictors, predictors, drop = FALSE]
  rxy <- R[predictors, outcome]

  ev <- eigen(Rxx, symmetric = TRUE)
  d <- ev$values
  if (min(d) <= tol * max(d)) {
    rw_stop("The predictors of '", label, "' are (almost) perfectly ",
            "collinear: the smallest eigenvalue of their correlation matrix ",
            "is ", format(min(d), digits = 3), ". Relative weights cannot be ",
            "computed.")
  }
  joint <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(joint) < -sqrt(.Machine$double.eps) * max(1, abs(joint))) {
    rw_stop("The correlation matrix of '", label, "' and its predictors is ",
            "not positive semi-definite (smallest eigenvalue ",
            format(min(joint), digits = 3), ").")
  }

  V <- ev$vectors
  Vr <- drop(crossprod(V, rxy))
  lambda <- V %*% (sqrt(d) * t(V))
  orth <- drop(V %*% (Vr / sqrt(d)))
  beta <- drop(V %*% (Vr / d))
  r2 <- sum(orth^2)
  if (r2 > 1 + 1e-8) {
    rw_stop("The R-squared of '", label, "' exceeds 1 (",
            format(r2, digits = 8), "); the correlation matrix is not valid.")
  }
  raw <- drop(lambda^2 %*% orth^2)
  if (r2 > 0) {
    rescaled <- 100 * raw / r2
  } else {
    rw_warn("The R-squared of '", label, "' is zero; rescaled weights are ",
            "undefined.")
    rescaled <- rep(NA_real_, length(raw))
  }

  names(raw) <- names(rescaled) <- names(orth) <- names(beta) <- predictors
  dimnames(lambda) <- list(predictors, predictors)
  list(raw = raw, rescaled = rescaled, r2 = r2, beta = beta, orth = orth,
       r = stats::setNames(unname(rxy), predictors), lambda = lambda, R = R)
}

rw_fit_specs <- function(specs, cors) {
  lapply(specs, function(s) {
    rw_compute(cors[[s$block]], s$outcome, s$predictors, label = s$label)
  })
}

# The statistics that are bootstrapped, in a fixed layout:
# for each regression, its raw weights, its rescaled weights and its R2.
rw_stat_vector <- function(fits) {
  unlist(lapply(fits, function(f) c(f$raw, f$rescaled, f$r2)),
         use.names = FALSE)
}

rw_stat_index <- function(specs) {
  out <- lapply(seq_along(specs), function(i) {
    p <- specs[[i]]$predictors
    data.frame(
      regression = i,
      statistic = c(rep("raw", length(p)), rep("rescaled", length(p)), "r2"),
      predictor = c(p, p, NA_character_),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

rw_sign <- function(x) as.integer(sign(x))

rw_percentile <- function(x, conf_level) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(NA_real_, NA_real_))
  }
  alpha <- (1 - conf_level) / 2
  stats::quantile(x, probs = c(alpha, 1 - alpha), type = 6, names = FALSE)
}

rw_label <- function(outcome, group = NA_character_, level = NA_character_) {
  n <- length(outcome)
  extra <- cbind(rep_len(as.character(group), n),
                 rep_len(as.character(level), n))
  vapply(seq_len(n), function(i) {
    e <- extra[i, ]
    e <- e[!is.na(e)]
    if (length(e)) {
      paste0(outcome[i], " [", paste(e, collapse = ", "), "]")
    } else {
      outcome[i]
    }
  }, character(1))
}

rw_collapse <- function(x) paste0("'", x, "'", collapse = ", ")

# Assemble the result object from specifications and fits.
rw_build <- function(specs, fits, info) {
  weights <- lapply(seq_along(specs), function(i) {
    s <- specs[[i]]
    f <- fits[[i]]
    data.frame(
      block = s$block, group = s$group, level = s$level, outcome = s$outcome,
      predictor = s$predictors, term = s$terms,
      raw = unname(f$raw), rescaled = unname(f$rescaled),
      beta = unname(f$beta), orth = unname(f$orth), r = unname(f$r),
      sign_beta = rw_sign(f$beta), sign_orth = rw_sign(f$orth),
      stringsAsFactors = FALSE, row.names = NULL
    )
  })
  r2 <- lapply(seq_along(specs), function(i) {
    s <- specs[[i]]
    data.frame(
      block = s$block, group = s$group, level = s$level, outcome = s$outcome,
      n = s$n, r2 = fits[[i]]$r2, r2_model = s$r2_model,
      stringsAsFactors = FALSE, row.names = NULL
    )
  })
  matrices <- lapply(fits, function(f) list(R = f$R, lambda = f$lambda))
  names(matrices) <- vapply(specs, function(s) s$label, character(1))

  structure(
    list(
      weights = do.call(rbind, weights),
      r2 = do.call(rbind, r2),
      matrices = matrices,
      boot = NULL,
      info = info
    ),
    class = "relweights"
  )
}

# Add percentile intervals from a replicate matrix `t` (one row per
# replicate, columns laid out as in rw_stat_index()).
rw_add_intervals <- function(x, t, t0, index, conf_level, n_nonadmissible = 0L,
                             method) {
  ok <- stats::complete.cases(t)
  n_failed <- sum(!ok)
  if (!any(ok)) {
    rw_stop("All ", nrow(t), " bootstrap replicates failed; no confidence ",
            "intervals can be computed.")
  }
  if (n_failed > 0L) {
    rw_warn(n_failed, " of ", nrow(t), " bootstrap replicates could not be ",
            "analysed (for example, a model did not converge or a ",
            "correlation matrix was singular or not positive ",
            "semi-definite) and were excluded.")
  }
  if (n_nonadmissible > 0L) {
    rw_warn(n_nonadmissible, " bootstrap replicates produced non-admissible ",
            "solutions; they were kept.")
  }
  ci <- vapply(seq_len(ncol(t)), function(j) {
    rw_percentile(t[ok, j], conf_level)
  }, numeric(2))

  is_raw <- index$statistic == "raw"
  is_rescaled <- index$statistic == "rescaled"
  is_r2 <- index$statistic == "r2"
  x$weights$raw_lower <- ci[1, is_raw]
  x$weights$raw_upper <- ci[2, is_raw]
  x$weights$rescaled_lower <- ci[1, is_rescaled]
  x$weights$rescaled_upper <- ci[2, is_rescaled]
  x$r2$r2_lower <- ci[1, is_r2]
  x$r2$r2_upper <- ci[2, is_r2]

  dimnames(t) <- NULL
  x$boot <- list(
    t = t, t0 = as.numeric(t0), index = index, n_boot = nrow(t),
    n_failed = n_failed,
    n_nonadmissible = as.integer(n_nonadmissible), conf_level = conf_level,
    method = method
  )
  x
}

rw_check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    rw_stop("`", name, "` must be TRUE or FALSE.")
  }
  invisible(x)
}

rw_check_boot_args <- function(n_boot, conf_level) {
  if (!is.numeric(n_boot) || length(n_boot) != 1L || !is.finite(n_boot) ||
      n_boot < 1 || n_boot != round(n_boot)) {
    rw_stop("`n_boot` must be a positive whole number.")
  }
  rw_check_conf_level(conf_level)
}

rw_check_outcome <- function(outcome) {
  if (!is.null(outcome) &&
      (!is.character(outcome) || !length(outcome) || anyNA(outcome))) {
    rw_stop("`outcome` must be a character vector of outcome names.")
  }
  invisible(outcome)
}

rw_warn_multi_column <- function(terms) {
  if (length(terms)) {
    rw_warn("Term(s) ", rw_collapse(terms), " produce several model-matrix ",
            "columns, each of which receives its own relative weight. For ",
            "factors, polynomials and splines these weights, and the ",
            "weights of the other predictors, depend on the contrasts or ",
            "basis used to code the term.")
  }
  invisible(terms)
}

rw_check_dots <- function(dots, allowed, what) {
  if (!length(dots)) {
    return(invisible(dots))
  }
  nms <- names(dots)
  if (is.null(nms) || any(!nzchar(nms))) {
    rw_stop("Arguments passed through `...` must be named.")
  }
  bad <- setdiff(nms, allowed)
  if (length(bad)) {
    rw_stop("Unknown argument(s) for relweights() with ", what, ": ",
            rw_collapse(bad), ". Arguments allowed in `...`: ",
            rw_collapse(allowed), ".")
  }
  invisible(dots)
}
