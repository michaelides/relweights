#' Print and extract relative weights
#'
#' `print()` shows, for every regression, \eqn{R^2}, the number of
#' observations and a table of weights, standardized coefficients,
#' correlations and signs. `as.data.frame()` returns one of the tables of
#' results.
#'
#' In the printed table, `raw` holds the raw weights, `rescaled` the
#' percentages of \eqn{R^2}, `beta` the standardized coefficients, `r` the
#' zero-order correlations, and `sign` the signs of `beta` and of the
#' orthogonal-counterpart coefficient `orth`, in that order (see
#' [relweights()]). Intervals are bootstrap confidence intervals (`CI`) or,
#' for Bayesian models, credible intervals (`CrI`). For Bayesian models, the
#' estimates are posterior means, and a second table gives the posterior
#' probabilities that `beta` and `orth` are positive and that each predictor
#' has the largest weight of its regression.
#'
#' @param x A `relweights` object.
#' @param digits Number of decimal places for raw weights, coefficients and
#'   correlations. Rescaled weights are shown with one decimal place.
#' @param ... Unused.
#'
#' @return `print()` returns `x` invisibly. `as.data.frame()` returns a data
#'   frame:
#' \describe{
#'   \item{`"weights"`, `"r2"`}{The `weights` or `r2` component of `x` (see
#'     the Value section of [relweights()]).}
#'   \item{`"ranks"`}{Bayesian models only. One row per predictor and rank,
#'     with the posterior probability (`probability`) that the predictor's
#'     raw weight has that rank in its regression (rank 1 is the largest
#'     weight).}
#'   \item{`"pairwise"`}{Bayesian models only. One row per ordered pair of
#'     predictors of a regression, with the posterior probability that the
#'     weight of `predictor` exceeds the weight of `other`.}
#' }
#'
#' @examples
#' rw <- relweights(lm(mpg ~ wt + hp + qsec, data = mtcars))
#' print(rw, digits = 2)
#' as.data.frame(rw)
#' as.data.frame(rw, what = "r2")
#' @export
print.relweights <- function(x, digits = 3, ...) {
  cat("Relative weight analysis (Johnson, 2000)\n")
  cat(strwrap(rw_model_description(x$info), width = getOption("width"),
              exdent = 2), sep = "\n")
  iv <- rw_interval_info(x)
  has_ci <- !is.null(iv)
  bayes <- !is.null(x$draws)
  w <- x$weights
  r2 <- x$r2
  sign_symbol <- function(s) c("-", "0", "+")[s + 2L]
  r2_model_label <- if (is.null(x$info$r2_model_label)) {
    "model"
  } else {
    x$info$r2_model_label
  }
  first <- if (bayes) rw_first_rank(x)

  for (i in seq_len(nrow(r2))) {
    rows <- which(w$block == r2$block[i] & w$outcome == r2$outcome[i])
    wi <- w[rows, , drop = FALSE]
    cat("\nOutcome: ", rw_label(r2$outcome[i], r2$group[i], r2$level[i]),
        "\n", sep = "")
    r2_text <- rw_fmt(r2$r2[i], digits)
    if (has_ci) {
      r2_text <- paste(r2_text, rw_fmt_ci(r2$r2_lower[i], r2$r2_upper[i],
                                          digits))
    }
    if (!is.na(r2$r2_model[i]) && abs(r2$r2[i] - r2$r2_model[i]) > 1e-6) {
      r2_text <- paste0(r2_text, " (", r2_model_label, ": ",
                        rw_fmt(r2$r2_model[i], digits), ")")
    }
    cat("R-squared = ", r2_text, ", n = ", r2$n[i], "\n\n", sep = "")

    tab <- data.frame(predictor = wi$predictor, raw = rw_fmt(wi$raw, digits),
                      check.names = FALSE, stringsAsFactors = FALSE)
    if (has_ci) {
      tab[[paste("raw", iv$label)]] <- rw_fmt_ci(wi$raw_lower, wi$raw_upper,
                                                 digits)
    }
    tab[["rescaled"]] <- rw_fmt(wi$rescaled, 1L)
    if (has_ci) {
      tab[[paste("rescaled", iv$label)]] <- rw_fmt_ci(wi$rescaled_lower,
                                                      wi$rescaled_upper, 1L)
    }
    tab[["beta"]] <- rw_fmt(wi$beta, digits)
    tab[["r"]] <- rw_fmt(wi$r, digits)
    tab[["sign"]] <- paste0(sign_symbol(wi$sign_beta), "/",
                            sign_symbol(wi$sign_orth))
    print(tab, row.names = FALSE, right = TRUE)

    if (bayes) {
      cat("\n")
      probs <- data.frame(predictor = wi$predictor,
                          `P(beta > 0)` = rw_fmt(wi$p_beta_pos, digits),
                          `P(orth > 0)` = rw_fmt(wi$p_orth_pos, digits),
                          `P(rank 1)` = rw_fmt(first[rows], digits),
                          check.names = FALSE, stringsAsFactors = FALSE)
      print(probs, row.names = FALSE, right = TRUE)
    }
  }

  cat("\nraw: raw weights (sum to R-squared); rescaled: % of R-squared;\n",
      "beta: standardized coefficients; r: zero-order correlations;\n",
      "sign: signs of beta / orthogonal-counterpart coefficients.\n", sep = "")
  differ <- w$sign_beta != w$sign_orth
  if (any(differ)) {
    where <- paste0(w$predictor[differ], " (",
                    rw_label(w$outcome[differ], w$group[differ],
                             w$level[differ]), ")")
    cat(strwrap(paste0("The two signs differ for: ",
                       paste(where, collapse = ", "), "."),
                width = getOption("width")), sep = "\n")
  }
  if (bayes) {
    text <- paste0("Estimates are posterior means. Intervals: ", iv$text,
                   " from ", iv$n, " posterior draws (", iv$nchains,
                   if (iv$nchains == 1L) " chain" else " chains", ")")
    if (x$draws$n_failed > 0L) {
      text <- paste0(text, "; ", x$draws$n_failed, " draws could not be ",
                     "analysed and were excluded")
    }
    cat(strwrap(paste0(text, "."), width = getOption("width")), sep = "\n")
    rhat <- w$rhat[is.finite(w$rhat)]
    ess <- w$ess_bulk[is.finite(w$ess_bulk)]
    if (length(rhat) && length(ess)) {
      text <- paste0("Raw weights: largest R-hat = ",
                     formatC(max(rhat), format = "f", digits = 3),
                     ", smallest bulk ESS = ", round(min(ess)), ".")
      if (max(rhat) > 1.01 || min(ess) < 400) {
        text <- paste(text, "Run longer chains for more reliable",
                      "results (R-hat should be below 1.01 and the ESS at",
                      "least 400).")
      }
      cat(strwrap(text, width = getOption("width")), sep = "\n")
    }
  } else if (has_ci) {
    b <- x$boot
    text <- paste0("Intervals: ", iv$text, " from ", b$n_boot, " replicates")
    if (b$n_failed > 0L) {
      text <- paste0(text, " (", b$n_failed, " failed and were excluded)")
    }
    cat(strwrap(paste0(text, "."), width = getOption("width")), sep = "\n")
  }
  invisible(x)
}

#' @rdname print.relweights
#' @param row.names,optional Standard arguments of [as.data.frame()];
#'   `row.names`, when not `NULL`, replaces the row names and `optional` is
#'   ignored.
#' @param what Which table to return: `"weights"` (one row per predictor),
#'   `"r2"` (one row per regression), or, for Bayesian models, `"ranks"` or
#'   `"pairwise"` (see Value).
#' @export
as.data.frame.relweights <- function(x, row.names = NULL, optional = FALSE,
                                     what = c("weights", "r2", "ranks",
                                              "pairwise"), ...) {
  what <- match.arg(what)
  out <- switch(what,
    weights = x$weights,
    r2 = x$r2,
    ranks = rw_rank_table(x),
    pairwise = rw_pairwise_table(x)
  )
  if (!is.null(row.names)) {
    rownames(out) <- row.names
  }
  out
}

rw_model_description <- function(info) {
  structure_text <- function(text) {
    if (isTRUE(info$ngroups > 1L)) {
      text <- paste0(text, "; groups: ",
                     paste(info$group_labels, collapse = ", "))
    }
    if (isTRUE(info$nlevels > 1L)) {
      text <- paste0(text, "; levels: ",
                     paste(info$level_labels, collapse = ", "))
    }
    text
  }
  text <- switch(
    if (is.null(info$model_class)) "lm" else info$model_class,
    lavaan = structure_text("lavaan model"),
    blavaan = structure_text("blavaan model"),
    brmsfit = paste0("brms model (gaussian, identity link): ",
                     paste(info$formula, collapse = "; ")),
    {
      text <- if (is.null(info$call)) {
        "lm"
      } else {
        paste(deparse(info$call, width.cutoff = 500L), collapse = " ")
      }
      if (isTRUE(info$weighted)) {
        text <- paste0(text, ", with prior weights")
      }
      text
    }
  )
  paste0("Model: ", text)
}

rw_fmt <- function(x, digits) {
  out <- formatC(x, format = "f", digits = digits)
  out[is.na(x)] <- "NA"
  out
}

rw_fmt_ci <- function(lower, upper, digits) {
  paste0("[", rw_fmt(lower, digits), ", ", rw_fmt(upper, digits), "]")
}
