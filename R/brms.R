#' @rdname relweights
#' @export
relweights.brmsfit <- function(object, outcome = NULL, ndraws = NULL,
                               conf_level = 0.95, ...) {
  rw_need("brms", "brms models")
  rw_need("posterior", "Bayesian models")
  rw_check_outcome(outcome)
  rw_check_bayes_args(ndraws, conf_level, list(...))

  parts <- rw_brms_parts(object)
  if (!is.null(outcome)) {
    available <- vapply(parts, function(p) p$label, character(1))
    resp <- vapply(parts, function(p) p$resp, character(1))
    unknown <- setdiff(outcome, c(available, resp))
    if (length(unknown)) {
      rw_stop("Outcome(s) ", rw_collapse(unknown), " are not responses of ",
              "the model. Available outcomes: ", rw_collapse(available), ".")
    }
    parts <- parts[available %in% outcome | resp %in% outcome]
  }

  sdata <- brms::standata(object)
  data <- lapply(parts, rw_brms_data, sdata = sdata)
  rw_warn_multi_column(unique(unlist(lapply(data, function(d) d$multi))))
  n_draws <- brms::ndraws(object)
  if (!length(n_draws) || n_draws < 1L) {
    rw_stop("The brms model contains no posterior draws.")
  }
  variables <- unlist(lapply(seq_along(parts), function(i) {
    c(paste0(parts[[i]]$b_prefix, data[[i]]$predictors), parts[[i]]$sigma)
  }))
  missing_vars <- setdiff(variables, posterior::variables(object))
  if (length(missing_vars)) {
    rw_stop("The posterior draws of ", rw_collapse(missing_vars), " were not ",
            "found in the brms model.")
  }
  arr <- unclass(posterior::as_draws_array(object, variable = variables))
  n_iter <- dim(arr)[1L]
  n_chains <- dim(arr)[2L]
  idx <- rw_thin_iterations(n_iter, n_chains, ndraws)
  flat <- function(v) as.vector(arr[idx, , v, drop = FALSE])
  chain <- rep(seq_len(n_chains), each = length(idx))
  iteration <- rep(idx, n_chains)

  specs <- vector("list", length(parts))
  draws <- vector("list", length(parts))
  for (i in seq_along(parts)) {
    p <- parts[[i]]
    d <- data[[i]]
    B <- vapply(paste0(p$b_prefix, d$predictors), flat,
                numeric(length(chain)))
    B <- matrix(B, ncol = length(d$predictors))
    draws[[i]] <- rw_fixed_x_draws(B, flat(p$sigma), d$X, d$y, p$label)
    specs[[i]] <- list(
      block = 1L, group = NA_character_, level = NA_character_,
      outcome = p$label, predictors = d$predictors, terms = d$terms,
      n = nrow(d$X), r2_model = NA_real_, label = p$label
    )
  }

  info <- list(
    model_class = "brmsfit", call = NULL, weighted = FALSE,
    ngroups = 1L, nlevels = 1L, group_labels = character(),
    level_labels = character(),
    multi_column_terms = unique(unlist(lapply(data, function(d) d$multi))),
    formula = vapply(parts, function(p) p$formula, character(1)),
    r2_model_label = "bayes_R2"
  )
  rw_bayes_build(specs, draws, chain, iteration, conf_level, info)
}

# Design matrix (without the intercept), response, predictor names and term
# labels of one response.
rw_brms_data <- function(p, sdata) {
  X <- sdata[[p$X]]
  y <- as.numeric(sdata[[p$Y]])
  if (is.null(X) || !length(y)) {
    rw_stop("The data of the model for '", p$label, "' could not be ",
            "retrieved with brms::standata().")
  }
  cols <- colnames(X)
  if (!"Intercept" %in% cols) {
    rw_stop("The model for '", p$label, "' has no intercept. Relative ",
            "weights are computed from correlations, which assume a model ",
            "with an intercept.")
  }
  keep <- cols != "Intercept"
  if (!any(keep)) {
    rw_stop("The model for '", p$label, "' has no predictors.")
  }
  assign <- attr(X, "assign")
  X <- unclass(X)[, keep, drop = FALSE]
  attr(X, "assign") <- NULL
  attr(X, "contrasts") <- NULL
  predictors <- colnames(X)
  terms <- predictors
  if (!is.null(assign) && length(assign) == length(cols)) {
    a <- assign[keep]
    known <- a >= 1L & a <= length(p$term_labels)
    terms[known] <- p$term_labels[a[known]]
  }
  counts <- table(terms)
  constant <- rw_constant_cols(X)
  if (length(constant)) {
    rw_stop("Predictor(s) ", rw_collapse(constant), " of '", p$label,
            "' have zero variance in the data used by the model.")
  }
  list(X = X, y = y, predictors = predictors, terms = terms,
       multi = names(counts)[counts > 1L])
}

# One element per response: its label, names of its draws and data, and the
# labels of its terms. Unsupported model features are errors.
rw_brms_parts <- function(object) {
  bterms <- brms::brmsterms(object$formula)
  multivariate <- inherits(bterms, "mvbrmsterms")
  parts <- if (multivariate) bterms$terms else list(bterms)
  lapply(parts, function(b) {
    label <- paste(deparse(b$respform[[2L]], width.cutoff = 500L),
                   collapse = " ")
    rw_brms_validate(b, label)
    resp <- if (multivariate) b$resp else ""
    suffix <- if (nzchar(resp)) paste0("_", resp) else ""
    fe <- b$dpars$mu$fe
    term_labels <- if (is.null(fe)) character() else {
      attr(stats::terms(fe), "term.labels")
    }
    rhs <- if (is.null(fe)) "1" else {
      paste(deparse(fe[[length(fe)]], width.cutoff = 500L), collapse = " ")
    }
    list(
      label = label, resp = resp,
      b_prefix = paste0("b", suffix, "_"), sigma = paste0("sigma", suffix),
      X = paste0("X", suffix), Y = paste0("Y", suffix),
      term_labels = term_labels,
      formula = paste(label, "~", sub("^1 \\+ ", "", rhs))
    )
  })
}

rw_brms_validate <- function(b, label) {
  family <- b$family$family
  link <- b$family$link
  if (!identical(family, "gaussian") || !identical(link, "identity")) {
    rw_stop("relweights() supports brms models with the gaussian family and ",
            "identity link; the model for '", label, "' uses the ", family,
            " family with the ", link, " link.")
  }
  if (length(b$adforms)) {
    rw_stop("The model for '", label, "' has addition terms (",
            paste0(names(b$adforms), "()", collapse = ", "), " after `|`), ",
            "which relative weights cannot take into account.")
  }
  dpars <- setdiff(names(b$dpars), "mu")
  if (length(dpars)) {
    rw_stop("The model for '", label, "' predicts the distributional ",
            "parameter(s) ", rw_collapse(dpars), ". Relative weights need a ",
            "constant residual variance; only the mean may depend on ",
            "predictors.")
  }
  if (length(b$nlpars) || inherits(b$dpars$mu, "btnl")) {
    rw_stop("The model for '", label, "' is non-linear. Relative weights ",
            "need a linear model.")
  }
  unsupported <- c(
    re = "group-level terms such as (1 | g)",
    sm = "smooth terms such as s() or t2()",
    sp = "special terms such as mo(), me() or mi()",
    gp = "Gaussian process terms (gp())",
    cs = "category-specific terms (cs())",
    offset = "an offset",
    ac = "autocorrelation terms such as ar() or arma()"
  )
  present <- vapply(names(unsupported), function(el) {
    z <- b$dpars$mu[[el]]
    !is.null(z) && NROW(z) > 0L
  }, logical(1))
  if (any(present)) {
    rw_stop("The model for '", label, "' contains ",
            paste(unsupported[present], collapse = " and "), ". Relative ",
            "weights support population-level effects only.")
  }
  invisible(TRUE)
}
