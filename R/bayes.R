# Shared machinery for Bayesian models (brms and blavaan).
#
# Each method computes, for every retained posterior draw and every
# regression, the raw and rescaled weights, standardized coefficients,
# orthogonal-counterpart coefficients, correlations, R2 and the model's own
# R2. rw_bayes_build() turns these per-draw results into the usual result
# object: posterior means in the main columns, equal-tailed credible
# intervals in the interval columns, and the draws in `x$draws`.

rw_need <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    rw_stop("Package '", pkg, "' is needed to analyse ", what, ". Install ",
            "it with install.packages(\"", pkg, "\").")
  }
  invisible(TRUE)
}

rw_check_conf_level <- function(conf_level) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    rw_stop("`conf_level` must be a single number between 0 and 1.")
  }
  invisible(TRUE)
}

rw_check_bayes_args <- function(ndraws, conf_level, dots) {
  if (length(dots)) {
    nms <- names(dots)
    if (is.null(nms)) nms <- rep("", length(dots))
    if (any(nms %in% c("bootstrap", "n_boot"))) {
      rw_stop("Bootstrap intervals are not used for Bayesian models: the ",
              "credible intervals come from the posterior draws. Use ",
              "`ndraws` and `conf_level` to control them.")
    }
    nms <- ifelse(nzchar(nms), nms, "<unnamed>")
    rw_stop("Unknown argument(s) for relweights() with Bayesian models: ",
            rw_collapse(nms), ".")
  }
  if (!is.null(ndraws) &&
      (!is.numeric(ndraws) || length(ndraws) != 1L || !is.finite(ndraws) ||
       ndraws < 1 || ndraws != round(ndraws))) {
    rw_stop("`ndraws` must be NULL or a positive whole number.")
  }
  rw_check_conf_level(conf_level)
}

# Iterations kept in every chain: all of them, or about ndraws / nchains
# evenly spaced ones (at least one per chain).
rw_thin_iterations <- function(n_iter, n_chains, ndraws = NULL) {
  if (is.null(ndraws) || ndraws >= n_iter * n_chains) {
    return(seq_len(n_iter))
  }
  k <- max(1L, floor(ndraws / n_chains))
  unique(round(seq(1, n_iter, length.out = k)))
}

# Per-draw weights for a regression whose predictors are fixed (conditioned
# on), as in brms: X is the design matrix without the intercept, B holds the
# coefficient draws (one row per draw) and sigma the residual SD draws.
#
# The predictor covariance S is the same in every draw, so its
# eigendecomposition (and Lambda) is computed once. In each draw the model
# implies var(y) = b'Sb + sigma^2, the standardized coefficients
# beta = b * sd(x) / sd(y) and r_xy = Rxx beta; since Lambda^2 = Rxx, the
# orthogonal-counterpart coefficients are Lambda^-1 r_xy = Lambda beta.
rw_fixed_x_draws <- function(B, sigma, X, y, outcome = "y", tol = 1e-10) {
  B <- as.matrix(B)
  X <- as.matrix(X)
  predictors <- colnames(X)
  S <- stats::cov(X)
  sdx <- sqrt(diag(S))
  Rxx <- stats::cov2cor(S)
  ev <- eigen(Rxx, symmetric = TRUE)
  d <- ev$values
  if (min(d) <= tol * max(d)) {
    rw_stop("The predictors of '", outcome, "' are (almost) perfectly ",
            "collinear: the smallest eigenvalue of their correlation matrix ",
            "is ", format(min(d), digits = 3), ". Relative weights cannot be ",
            "computed.")
  }
  lambda <- ev$vectors %*% (sqrt(d) * t(ev$vectors))
  dimnames(lambda) <- list(predictors, predictors)

  var_fit <- rowSums((B %*% S) * B)
  var_y <- var_fit + sigma^2
  beta <- sweep(B, 2L, sdx, `*`) / sqrt(var_y)
  orth <- beta %*% lambda
  raw <- orth^2 %*% lambda^2
  r2 <- rowSums(orth^2)
  rescaled <- 100 * raw / r2
  r <- beta %*% Rxx

  # The residual-based R2 of Gelman et al. (2019), as in brms::bayes_R2():
  # var(fitted) / (var(fitted) + var(y - fitted)) within each draw.
  cxy <- drop(stats::cov(X, y))
  var_e <- stats::var(y) - 2 * drop(B %*% cxy) + var_fit
  r2_model <- var_fit / (var_fit + var_e)

  out <- list(raw = raw, rescaled = rescaled, beta = beta, orth = orth,
              r = r)
  out <- lapply(out, function(m) {
    m <- unname(as.matrix(m))
    colnames(m) <- predictors
    m
  })
  r_mean <- colMeans(out$r)
  R <- rbind(c(1, r_mean), cbind(r_mean, Rxx))
  dimnames(R) <- list(c(outcome, predictors), c(outcome, predictors))
  c(out, list(r2 = unname(r2), r2_model = unname(r2_model), R = R,
              lambda = lambda))
}

rw_quantiles <- function(v, conf_level) {
  alpha <- (1 - conf_level) / 2
  stats::quantile(v, probs = c(alpha, 1 - alpha), names = FALSE, type = 7)
}

# R-hat and bulk ESS of one quantity, from draws stored chain by chain.
rw_convergence <- function(v, chain) {
  na <- c(rhat = NA_real_, ess_bulk = NA_real_)
  counts <- tabulate(match(chain, unique(chain)))
  if (anyNA(v) || length(unique(counts)) != 1L || counts[1L] < 4L) {
    return(na)
  }
  m <- matrix(v, ncol = length(counts))
  out <- tryCatch(
    suppressWarnings(c(rhat = posterior::rhat(m),
                       ess_bulk = posterior::ess_bulk(m))),
    error = function(e) na
  )
  out[!is.finite(out)] <- NA_real_
  out
}

# Build the result object from per-draw results.
#
# `draws` has one element per regression (in the order of `specs`), each a
# list with matrices raw, rescaled, beta, orth and r (draws x predictors),
# vectors r2 and r2_model, and the matrices R and lambda used for
# `x$matrices`. Rows with missing values are failed draws.
rw_bayes_build <- function(specs, draws, chain, iteration, conf_level, info) {
  ok <- Reduce(`&`, lapply(draws, function(d) {
    stats::complete.cases(d$raw, d$rescaled, d$r2)
  }))
  n_failed <- sum(!ok)
  if (!any(ok)) {
    rw_stop("None of the ", length(ok), " posterior draws could be ",
            "analysed.")
  }
  if (n_failed > 0L) {
    rw_warn(n_failed, " of ", length(ok), " posterior draws could not be ",
            "analysed (for example, a correlation matrix was singular or ",
            "not positive semi-definite) and were excluded.")
  }

  keep <- function(m) m[ok, , drop = FALSE]
  fits <- lapply(draws, function(d) {
    list(raw = colMeans(keep(d$raw)), rescaled = colMeans(keep(d$rescaled)),
         beta = colMeans(keep(d$beta)), orth = colMeans(keep(d$orth)),
         r = colMeans(keep(d$r)), r2 = mean(d$r2[ok]), R = d$R,
         lambda = d$lambda)
  })
  for (i in seq_along(specs)) {
    specs[[i]]$r2_model <- mean(draws[[i]]$r2_model[ok])
  }
  info$intervals <- "credible"
  info$ndraws <- sum(ok)
  info$nchains <- length(unique(chain))
  x <- rw_build(specs, fits, info)

  interval <- function(name) {
    do.call(cbind, lapply(draws, function(d) {
      apply(keep(d[[name]]), 2L, rw_quantiles, conf_level = conf_level)
    }))
  }
  ci_raw <- interval("raw")
  ci_rescaled <- interval("rescaled")
  x$weights$raw_lower <- ci_raw[1L, ]
  x$weights$raw_upper <- ci_raw[2L, ]
  x$weights$rescaled_lower <- ci_rescaled[1L, ]
  x$weights$rescaled_upper <- ci_rescaled[2L, ]
  ci_r2 <- vapply(draws, function(d) rw_quantiles(d$r2[ok], conf_level),
                  numeric(2))
  x$r2$r2_lower <- ci_r2[1L, ]
  x$r2$r2_upper <- ci_r2[2L, ]

  x$weights$p_beta_pos <- unlist(lapply(draws, function(d) {
    colMeans(keep(d$beta) > 0)
  }), use.names = FALSE)
  x$weights$p_orth_pos <- unlist(lapply(draws, function(d) {
    colMeans(keep(d$orth) > 0)
  }), use.names = FALSE)
  conv <- do.call(cbind, lapply(draws, function(d) {
    vapply(seq_len(ncol(d$raw)), function(j) rw_convergence(d$raw[, j], chain),
           numeric(2))
  }))
  x$weights$rhat <- conv[1L, ]
  x$weights$ess_bulk <- conv[2L, ]
  if (any(x$weights$rhat > 1.05, na.rm = TRUE)) {
    rw_warn("The posterior draws of some relative weights have not ",
            "converged (largest R-hat ",
            format(max(x$weights$rhat, na.rm = TRUE), digits = 3), "). Run ",
            "longer chains before interpreting the results.")
  }

  labels <- vapply(specs, function(s) s$label, character(1))
  x$draws <- list(
    chain = chain, iteration = iteration, ok = ok, conf_level = conf_level,
    n_failed = n_failed,
    regressions = stats::setNames(lapply(draws, function(d) {
      d[c("raw", "rescaled", "beta", "orth", "r", "r2", "r2_model")]
    }), labels)
  )
  x
}

# Draws of the raw weights of one regression (rows = successful draws).
rw_regression_draws <- function(x, i) {
  label <- rw_label(x$r2$outcome[i], x$r2$group[i], x$r2$level[i])
  x$draws$regressions[[label]]$raw[x$draws$ok, , drop = FALSE]
}

rw_require_draws <- function(x, what) {
  if (is.null(x$draws)) {
    rw_stop("`what = \"", what, "\"` needs the posterior draws of a ",
            "Bayesian model (brms or blavaan).")
  }
  invisible(TRUE)
}

# P(rank of each predictor = k), with rank 1 for the largest weight.
rw_rank_table <- function(x) {
  rw_require_draws(x, "ranks")
  out <- lapply(seq_len(nrow(x$r2)), function(i) {
    d <- rw_regression_draws(x, i)
    p <- ncol(d)
    ranks <- if (p == 1L) {
      matrix(1L, nrow(d), 1L)
    } else {
      t(apply(-d, 1L, rank, ties.method = "first"))
    }
    probs <- vapply(seq_len(p), function(j) {
      tabulate(ranks[, j], nbins = p) / nrow(d)
    }, numeric(p))
    data.frame(
      block = x$r2$block[i], group = x$r2$group[i], level = x$r2$level[i],
      outcome = x$r2$outcome[i],
      predictor = rep(colnames(d), each = p),
      rank = rep(seq_len(p), times = p),
      probability = as.vector(probs),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

# P(weight of predictor > weight of other) for every ordered pair.
rw_pairwise_table <- function(x) {
  rw_require_draws(x, "pairwise")
  out <- lapply(seq_len(nrow(x$r2)), function(i) {
    d <- rw_regression_draws(x, i)
    p <- ncol(d)
    if (p < 2L) {
      return(NULL)
    }
    pairs <- expand.grid(k = seq_len(p), j = seq_len(p))
    pairs <- pairs[pairs$j != pairs$k, , drop = FALSE]
    data.frame(
      block = x$r2$block[i], group = x$r2$group[i], level = x$r2$level[i],
      outcome = x$r2$outcome[i],
      predictor = colnames(d)[pairs$j], other = colnames(d)[pairs$k],
      probability = vapply(seq_len(nrow(pairs)), function(r) {
        mean(d[, pairs$j[r]] > d[, pairs$k[r]])
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  if (is.null(out)) {
    rw_stop("Pairwise probabilities need regressions with at least two ",
            "predictors.")
  }
  rownames(out) <- NULL
  out
}

# Probability that each predictor has the largest weight of its regression.
rw_first_rank <- function(x) {
  unlist(lapply(seq_len(nrow(x$r2)), function(i) {
    d <- rw_regression_draws(x, i)
    best <- max.col(d, ties.method = "first")
    tabulate(best, nbins = ncol(d)) / nrow(d)
  }), use.names = FALSE)
}

# Description of the intervals in a result, or NULL without intervals.
rw_interval_info <- function(x) {
  if (!is.null(x$boot)) {
    return(list(
      type = "bootstrap", label = "CI", conf_level = x$boot$conf_level,
      n = x$boot$n_boot - x$boot$n_failed,
      text = paste0(format(100 * x$boot$conf_level), "% percentile ",
                    "bootstrap intervals")
    ))
  }
  if (!is.null(x$draws)) {
    return(list(
      type = "credible", label = "CrI", conf_level = x$draws$conf_level,
      n = sum(x$draws$ok), nchains = length(unique(x$draws$chain)),
      text = paste0(format(100 * x$draws$conf_level), "% equal-tailed ",
                    "credible intervals")
    ))
  }
  NULL
}
