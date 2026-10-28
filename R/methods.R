#' Print and extract relative weights
#'
#' `print()` shows, for every regression, \eqn{R^2}, the number of
#' observations and a table of weights, standardized coefficients,
#' correlations and signs. `as.data.frame()` returns the table of weights or
#' the table of \eqn{R^2} values.
#'
#' In the printed table, `raw` holds the raw weights, `rescaled` the
#' percentages of \eqn{R^2}, `beta` the standardized coefficients, `r` the
#' zero-order correlations, and `sign` the signs of `beta` and of the
#' orthogonal-counterpart coefficient `orth`, in that order (see
#' [relweights()]).
#'
#' @param x A `relweights` object.
#' @param digits Number of decimal places for raw weights, coefficients and
#'   correlations. Rescaled weights are shown with one decimal place.
#' @param ... Unused.
#'
#' @return `print()` returns `x` invisibly. `as.data.frame()` returns the
#'   `weights` or the `r2` component of `x` (see the Value section of
#'   [relweights()]).
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
  has_ci <- !is.null(x$boot)
  w <- x$weights
  r2 <- x$r2
  sign_symbol <- function(s) c("-", "0", "+")[s + 2L]

  for (i in seq_len(nrow(r2))) {
    wi <- w[w$block == r2$block[i] & w$outcome == r2$outcome[i], ,
            drop = FALSE]
    cat("\nOutcome: ", rw_label(r2$outcome[i], r2$group[i], r2$level[i]),
        "\n", sep = "")
    r2_text <- rw_fmt(r2$r2[i], digits)
    if (has_ci) {
      r2_text <- paste(r2_text, rw_fmt_ci(r2$r2_lower[i], r2$r2_upper[i],
                                          digits))
    }
    if (!is.na(r2$r2_model[i]) && abs(r2$r2[i] - r2$r2_model[i]) > 1e-6) {
      r2_text <- paste0(r2_text, " (model: ", rw_fmt(r2$r2_model[i], digits),
                        ")")
    }
    cat("R-squared = ", r2_text, ", n = ", r2$n[i], "\n\n", sep = "")

    tab <- data.frame(predictor = wi$predictor, raw = rw_fmt(wi$raw, digits),
                      check.names = FALSE, stringsAsFactors = FALSE)
    if (has_ci) {
      tab[["raw CI"]] <- rw_fmt_ci(wi$raw_lower, wi$raw_upper, digits)
    }
    tab[["rescaled"]] <- rw_fmt(wi$rescaled, 1L)
    if (has_ci) {
      tab[["rescaled CI"]] <- rw_fmt_ci(wi$rescaled_lower, wi$rescaled_upper,
                                        1L)
    }
    tab[["beta"]] <- rw_fmt(wi$beta, digits)
    tab[["r"]] <- rw_fmt(wi$r, digits)
    tab[["sign"]] <- paste0(sign_symbol(wi$sign_beta), "/",
                            sign_symbol(wi$sign_orth))
    print(tab, row.names = FALSE, right = TRUE)
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
  if (has_ci) {
    b <- x$boot
    text <- paste0("Intervals: ", format(100 * b$conf_level), "% percentile ",
                   "bootstrap intervals from ", b$n_boot, " replicates")
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
#' @param what Which table to return: `"weights"` (one row per predictor) or
#'   `"r2"` (one row per regression).
#' @export
as.data.frame.relweights <- function(x, row.names = NULL, optional = FALSE,
                                     what = c("weights", "r2"), ...) {
  what <- match.arg(what)
  out <- x[[what]]
  if (!is.null(row.names)) {
    rownames(out) <- row.names
  }
  out
}

rw_model_description <- function(info) {
  if (identical(info$model_class, "lavaan")) {
    text <- "lavaan model"
    if (isTRUE(info$ngroups > 1L)) {
      text <- paste0(text, "; groups: ",
                     paste(info$group_labels, collapse = ", "))
    }
    if (isTRUE(info$nlevels > 1L)) {
      text <- paste0(text, "; levels: ",
                     paste(info$level_labels, collapse = ", "))
    }
  } else {
    text <- if (is.null(info$call)) {
      "lm"
    } else {
      paste(deparse(info$call, width.cutoff = 500L), collapse = " ")
    }
    if (isTRUE(info$weighted)) {
      text <- paste0(text, ", with prior weights")
    }
  }
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
