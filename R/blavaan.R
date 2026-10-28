#' @rdname relweights
#' @export
relweights.blavaan <- function(object, outcome = NULL, ndraws = NULL,
                               conf_level = 0.95, ...) {
  rw_need("blavaan", "blavaan models")
  rw_need("posterior", "Bayesian models")
  rw_check_outcome(outcome)
  rw_check_bayes_args(ndraws, conf_level, list(...))

  specs <- rw_lavaan_specs(object, outcome)
  rw_lavaan_check_vars(specs, rw_lavaan_cors(object))
  mcmc <- rw_blavaan_draws(object, ndraws)
  resvar <- rw_blavaan_resvar(object, specs)

  n_draws <- nrow(mcmc$x)
  draws <- lapply(specs, function(s) {
    p <- length(s$predictors)
    empty <- matrix(NA_real_, n_draws, p, dimnames = list(NULL, s$predictors))
    list(raw = empty, rescaled = empty, beta = empty, orth = empty, r = empty,
         r2 = rep(NA_real_, n_draws), r2_model = rep(NA_real_, n_draws),
         R = 0, lambda = 0)
  })
  n_ok <- 0L
  for (d in seq_len(n_draws)) {
    res <- tryCatch(
      suppressWarnings(rw_blavaan_draw(object, specs, mcmc$x[d, ], resvar)),
      error = function(e) NULL
    )
    if (is.null(res)) {
      next
    }
    n_ok <- n_ok + 1L
    for (i in seq_along(specs)) {
      f <- res$fits[[i]]
      draws[[i]]$raw[d, ] <- f$raw
      draws[[i]]$rescaled[d, ] <- f$rescaled
      draws[[i]]$beta[d, ] <- f$beta
      draws[[i]]$orth[d, ] <- f$orth
      draws[[i]]$r[d, ] <- f$r
      draws[[i]]$r2[d] <- f$r2
      draws[[i]]$r2_model[d] <- res$r2_model[i]
      draws[[i]]$R <- draws[[i]]$R + f$R
      draws[[i]]$lambda <- draws[[i]]$lambda + f$lambda
    }
  }
  if (n_ok > 0L) {
    for (i in seq_along(draws)) {
      draws[[i]]$R <- draws[[i]]$R / n_ok
      draws[[i]]$lambda <- draws[[i]]$lambda / n_ok
    }
  }

  info <- list(
    model_class = "blavaan", call = NULL, weighted = NA,
    ngroups = as.integer(lavaan::lavInspect(object, "ngroups")),
    nlevels = as.integer(lavaan::lavInspect(object, "nlevels")),
    group_labels = as.character(lavaan::lavInspect(object, "group.label")),
    level_labels = as.character(lavaan::lavInspect(object, "level.label")),
    multi_column_terms = character(),
    r2_model_label = "model"
  )
  x <- rw_bayes_build(specs, draws, mcmc$chain, mcmc$iteration, conf_level,
                      info)

  ok <- x$draws$ok
  differs <- vapply(draws, function(d) {
    any(is.finite(d$r2_model[ok])) &&
      mean(abs(d$r2[ok] - d$r2_model[ok]), na.rm = TRUE) > 1e-6
  }, logical(1))
  rw_warn_implied(specs[differs])
  x
}

# Posterior draws of the free parameters, thinned and stacked chain by chain.
rw_blavaan_draws <- function(object, ndraws) {
  chains <- tryCatch(blavaan::blavInspect(object, "mcmc"),
                     error = function(e) NULL)
  if (is.null(chains) || !length(chains)) {
    rw_stop("The blavaan model contains no posterior draws.")
  }
  chains <- lapply(chains, as.matrix)
  free <- names(lavaan::coef(object))
  n_free <- length(free)
  ok <- all(vapply(chains, function(m) {
    ncol(m) == n_free && identical(colnames(m), free)
  }, logical(1)))
  if (!ok) {
    rw_stop("The posterior draws of the blavaan model do not match its free ",
            "parameters (", n_free, " expected), so they cannot be assigned ",
            "to the model.")
  }
  n_iter <- vapply(chains, nrow, integer(1))
  if (length(unique(n_iter)) != 1L || n_iter[1L] < 1L) {
    rw_stop("The chains of the blavaan model have different lengths.")
  }
  idx <- rw_thin_iterations(n_iter[1L], length(chains), ndraws)
  list(
    x = do.call(rbind, lapply(chains, function(m) m[idx, , drop = FALSE])),
    chain = rep(seq_along(chains), each = length(idx)),
    iteration = rep(idx, length(chains))
  )
}

# Where to find the residual variance of each outcome: the position of its
# free parameter, or its fixed value.
rw_blavaan_resvar <- function(object, specs) {
  pt <- lavaan::parTable(object)
  lapply(specs, function(s) {
    row <- which(pt$op == "~~" & pt$lhs == s$outcome &
                   pt$rhs == s$outcome & pt$block == s$block)
    if (length(row) != 1L) {
      return(list(free = 0L, est = NA_real_))
    }
    list(free = as.integer(pt$free[row]), est = as.numeric(pt$est[row]))
  })
}

# Relative weights and the model's R2 for one parameter vector `x`.
rw_blavaan_draw <- function(object, specs, x, resvar) {
  fit <- object
  fit@Model <- lavaan::lav_model_set_parameters(object@Model, x = unname(x))
  covs <- lapply(rw_as_list(lavaan::lavInspect(fit, "cov.all")), unclass)
  cors <- lapply(covs, stats::cov2cor)
  fits <- rw_fit_specs(specs, cors)
  r2_model <- vapply(seq_along(specs), function(i) {
    s <- specs[[i]]
    rv <- resvar[[i]]
    value <- if (rv$free > 0L) x[[rv$free]] else rv$est
    1 - value / covs[[s$block]][s$outcome, s$outcome]
  }, numeric(1))
  list(fits = fits, r2_model = r2_model, cors = cors)
}
