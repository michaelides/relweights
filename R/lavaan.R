#' @rdname relweights
#' @export
relweights.lavaan <- function(object, outcome = NULL, bootstrap = FALSE,
                              n_boot = 1000, conf_level = 0.95, ...) {
  if (!identical(class(object)[1L], "lavaan")) {
    rw_stop("relweights() supports models fitted with lavaan; objects of ",
            "class ", rw_collapse(class(object)[1L]), " are not supported.")
  }
  if (!is.null(outcome) &&
      (!is.character(outcome) || !length(outcome) || anyNA(outcome))) {
    rw_stop("`outcome` must be a character vector of outcome names.")
  }
  rw_check_flag(bootstrap, "bootstrap")
  rw_check_boot_args(n_boot, conf_level)
  dots <- list(...)
  rw_check_dots(dots, c("parallel", "ncpus", "cl", "iseed", "type"),
                "lavaan models")

  if (!isTRUE(lavaan::lavInspect(object, "converged"))) {
    rw_stop("The lavaan model has not converged, so its implied ",
            "correlations cannot be trusted.")
  }
  admissible <- tryCatch(
    suppressWarnings(lavaan::lavInspect(object, "post.check")),
    error = function(e) NA
  )
  if (isFALSE(admissible)) {
    rw_warn("The lavaan solution is not admissible (for example, it has ",
            "negative variances); interpret the relative weights with ",
            "caution.")
  }

  specs <- rw_lavaan_specs(object, outcome)
  cors <- rw_lavaan_cors(object)
  rw_lavaan_check_vars(specs, cors)
  fits <- rw_fit_specs(specs, cors)
  rw_lavaan_crosscheck(object, specs, fits)

  call <- tryCatch(lavaan::lavInspect(object, "call"),
                   error = function(e) NULL)
  if (is.list(call)) {
    call <- as.call(call)
  }
  info <- list(
    model_class = "lavaan", call = call, weighted = NA,
    ngroups = as.integer(lavaan::lavInspect(object, "ngroups")),
    nlevels = as.integer(lavaan::lavInspect(object, "nlevels")),
    group_labels = as.character(lavaan::lavInspect(object, "group.label")),
    level_labels = as.character(lavaan::lavInspect(object, "level.label")),
    multi_column_terms = character()
  )
  out <- rw_build(specs, fits, info)
  if (bootstrap) {
    out <- rw_boot_lavaan(out, object, specs, fits, n_boot, conf_level, dots)
  }
  out
}

rw_as_list <- function(x) {
  if (is.list(x)) x else list(x)
}

rw_lavaan_cors <- function(object) {
  lapply(rw_as_list(lavaan::lavInspect(object, "cor.all")), unclass)
}

# One specification per (block, outcome) pair, in parameter-table order.
rw_lavaan_specs <- function(object, outcome = NULL) {
  pt <- lavaan::parTable(object)
  reg <- pt[pt$op == "~", c("lhs", "rhs", "block"), drop = FALSE]
  if (!nrow(reg)) {
    rw_stop("The lavaan model contains no regressions (`~`).")
  }
  if (!is.null(outcome)) {
    unknown <- setdiff(outcome, reg$lhs)
    if (length(unknown)) {
      rw_stop("Outcome(s) ", rw_collapse(unknown), " are not regressed on ",
              "any predictor in the model. Available outcomes: ",
              rw_collapse(unique(reg$lhs)), ".")
    }
    reg <- reg[reg$lhs %in% outcome, , drop = FALSE]
  }

  ngroups <- lavaan::lavInspect(object, "ngroups")
  nlevels <- lavaan::lavInspect(object, "nlevels")
  group_labels <- lavaan::lavInspect(object, "group.label")
  level_labels <- lavaan::lavInspect(object, "level.label")
  nobs <- unlist(lavaan::lavInspect(object, "nobs"))
  nclusters <- if (nlevels > 1L) {
    unlist(lavaan::lavInspect(object, "nclusters"))
  }
  rsquare <- rw_as_list(lavaan::lavInspect(object, "rsquare"))

  key <- paste(reg$block, reg$lhs, sep = "\r")
  lapply(unique(key), function(k) {
    rows <- reg[key == k, , drop = FALSE]
    b <- as.integer(rows$block[1L])
    y <- rows$lhs[1L]
    g <- (b - 1L) %/% nlevels + 1L
    l <- (b - 1L) %% nlevels + 1L
    group <- if (ngroups > 1L) group_labels[g] else NA_character_
    level <- if (nlevels > 1L) level_labels[l] else NA_character_
    n <- if (l == 1L) nobs[g] else nclusters[g]
    r2_model <- if (length(rsquare) >= b) unname(rsquare[[b]][y]) else NA
    predictors <- unique(rows$rhs)
    list(
      block = b, group = unname(group), level = unname(level), outcome = y,
      predictors = predictors, terms = predictors,
      n = as.integer(unname(n)), r2_model = as.numeric(r2_model),
      label = rw_label(y, unname(group), unname(level))
    )
  })
}

rw_lavaan_check_vars <- function(specs, cors) {
  for (s in specs) {
    absent <- setdiff(c(s$outcome, s$predictors), rownames(cors[[s$block]]))
    if (length(absent)) {
      rw_stop("Variable(s) ", rw_collapse(absent), " of the regression of '",
              s$label, "' are not part of the model-implied correlation ",
              "matrix. This happens when lavaan conditions on exogenous ",
              "covariates (conditional.x = TRUE, the default when endogenous ",
              "variables are ordered); refit the model with ",
              "conditional.x = FALSE.")
    }
  }
  invisible(TRUE)
}

# Compare the implied regression with the model's own regression.
rw_lavaan_crosscheck <- function(object, specs, fits, tol = 1e-6) {
  pt <- lavaan::parTable(object)
  covs <- rw_as_list(lavaan::lavInspect(object, "cov.all"))
  differs <- vapply(seq_along(specs), function(i) {
    s <- specs[[i]]
    f <- fits[[i]]
    sds <- sqrt(diag(unclass(covs[[s$block]])))
    rows <- pt$op == "~" & pt$block == s$block & pt$lhs == s$outcome
    est <- stats::setNames(pt$est[rows], pt$rhs[rows])[s$predictors]
    std <- est * sds[s$predictors] / sds[s$outcome]
    beta_differs <- any(!is.finite(std)) || any(abs(std - f$beta) > tol)
    r2_differs <- !is.na(s$r2_model) && abs(f$r2 - s$r2_model) > tol
    beta_differs || r2_differs
  }, logical(1))
  if (any(differs)) {
    labels <- vapply(specs[differs], function(s) s$label, character(1))
    rw_warn("For ", rw_collapse(labels), ", the regression implied by the ",
            "model-implied correlations differs from the model's own ",
            "regression, because the outcome's residual is correlated with ",
            "its predictors (for example through a residual covariance or a ",
            "feedback loop). The relative weights, `beta` and `r2` describe ",
            "the implied regression; `r2_model` holds the R-squared reported ",
            "by lavaan.")
  }
  invisible(differs)
}

rw_pick_formal <- function(formal_names, candidates) {
  hit <- candidates[candidates %in% formal_names]
  if (!length(hit)) {
    rw_stop("This version of lavaan::lavBootstrap() is not supported: none ",
            "of its arguments is named ", rw_collapse(candidates), ".")
  }
  hit[1L]
}

rw_boot_lavaan <- function(x, object, specs, fits, n_boot, conf_level, dots) {
  if (lavaan::lavInspect(object, "nlevels") > 1L) {
    rw_stop("Bootstrap confidence intervals are not available for ",
            "multilevel lavaan models.")
  }
  data_type <- tryCatch(object@Data@data.type, error = function(e) NA)
  if (identical(data_type, "moment") || identical(data_type, "none")) {
    rw_stop("Bootstrap confidence intervals need a lavaan model fitted to ",
            "raw data, not to summary statistics.")
  }

  t0 <- rw_stat_vector(fits)
  n_stat <- length(t0)
  fun <- function(fit) {
    tryCatch(
      rw_stat_vector(suppressWarnings(rw_fit_specs(specs,
                                                   rw_lavaan_cors(fit)))),
      error = function(e) rep(NA_real_, n_stat)
    )
  }

  formal_names <- names(formals(lavaan::lavBootstrap))
  args <- list(object)
  args[[rw_pick_formal(formal_names, c("R", "r"))]] <- n_boot
  args[[rw_pick_formal(formal_names, c("FUN", "fun"))]] <- fun
  args <- c(args, dots)

  reps <- tryCatch(
    withCallingHandlers(
      do.call(lavaan::lavBootstrap, args),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) {
      rw_stop("Bootstrapping the lavaan model failed: ", conditionMessage(e))
    }
  )
  n_nonadmissible <- length(attr(reps, "nonadmissible"))
  reps <- matrix(as.numeric(reps), nrow = nrow(reps), ncol = ncol(reps))
  rw_add_intervals(x, reps, t0, rw_stat_index(specs), conf_level,
                   n_nonadmissible = n_nonadmissible,
                   method = "lavaan::lavBootstrap")
}
