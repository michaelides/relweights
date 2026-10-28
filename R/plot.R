#' Plot relative weights
#'
#' `autoplot()` draws the results of [relweights()] with 'ggplot2' and
#' returns the plot; `plot()` draws it and returns it invisibly.
#'
#' Four plots are available:
#' \describe{
#'   \item{`"weights"`}{A bar chart of the relative weights of each
#'     regression, sorted from the largest to the smallest weight, with one
#'     panel per regression. With `sign = "beta"` or `sign = "orthogonal"`,
#'     each weight takes the sign of the standardized coefficient or of the
#'     orthogonal-counterpart coefficient (see [relweights()]). Bootstrap
#'     percentile intervals or, for Bayesian models, credible intervals are
#'     drawn as error bars when available.}
#'   \item{`"r2"`}{One bar per regression, divided into the raw weights of
#'     its predictors, so that the bar length is \eqn{R^2}. Bootstrap or
#'     credible intervals for \eqn{R^2} are drawn when available.}
#'   \item{`"compare"`}{The relative weights next to two simpler indices of
#'     importance, the squared standardized coefficients and the squared
#'     zero-order correlations, each expressed as a percentage of its total
#'     over the predictors of the regression.}
#'   \item{`"groups"`}{A dot plot comparing the weight of each predictor
#'     across the groups and levels of a 'lavaan' or 'blavaan' model, with
#'     one panel per outcome. Requires a model with more than one group or
#'     level.}
#' }
#'
#' For Bayesian models, bars and points show posterior means.
#'
#' @param object,x A `relweights` object.
#' @param type The plot to draw: `"weights"`, `"r2"`, `"compare"` or
#'   `"groups"` (see Details).
#' @param scale Plot `"rescaled"` weights (percentages of \eqn{R^2}) or
#'   `"raw"` weights. Used by the `"weights"` and `"groups"` plots.
#' @param sign `"none"`, `"beta"` or `"orthogonal"`: give each weight the
#'   sign of its standardized coefficient (`"beta"`) or of its
#'   orthogonal-counterpart coefficient (`"orthogonal"`). Used by the
#'   `"weights"` and `"groups"` plots.
#' @param ci Logical. Draw bootstrap or credible intervals when they are
#'   available?
#' @param outcome Optional character vector of outcomes to plot.
#' @param ... For `autoplot()`, unused (must be empty). For `plot()`,
#'   arguments passed to `autoplot()`.
#'
#' @return `autoplot()` returns a ggplot object; `plot()` returns the same
#'   object invisibly.
#'
#' @seealso [relweights()]
#'
#' @examples
#' fit <- lm(mpg ~ cyl + disp + hp + wt + qsec, data = mtcars)
#' rw <- relweights(fit)
#' autoplot(rw)
#' autoplot(rw, sign = "beta")
#' autoplot(rw, type = "r2")
#' autoplot(rw, type = "compare")
#'
#' hs_model <- "
#'   visual  =~ x1 + x2 + x3
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#'   speed ~ visual + textual + ageyr
#' "
#' hs_fit <- lavaan::sem(hs_model, data = lavaan::HolzingerSwineford1939,
#'                       group = "school")
#' autoplot(relweights(hs_fit), type = "groups")
#' @export
autoplot.relweights <- function(object,
                                type = c("weights", "r2", "compare", "groups"),
                                scale = c("rescaled", "raw"),
                                sign = c("none", "beta", "orthogonal"),
                                ci = TRUE, outcome = NULL, ...) {
  if (...length()) {
    rw_stop("Unused argument(s) in autoplot(): ",
            paste(rw_dot_names(...), collapse = ", "), ".")
  }
  type <- match.arg(type)
  scale <- match.arg(scale)
  sign <- match.arg(sign)
  rw_check_flag(ci, "ci")
  x <- rw_subset_outcome(object, outcome)
  switch(type,
    weights = rw_plot_weights(x, scale, sign, ci),
    r2 = rw_plot_r2(x, ci),
    compare = rw_plot_compare(x),
    groups = rw_plot_groups(x, scale, sign, ci)
  )
}

#' @rdname autoplot.relweights
#' @export
plot.relweights <- function(x, ...) {
  p <- autoplot.relweights(x, ...)
  print(p)
  invisible(p)
}

# Helpers -----------------------------------------------------------------

rw_palette <- c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00",
                "#CC79A7", "#F0E442", "#999999")

rw_dot_names <- function(...) {
  nms <- names(list(...))
  if (is.null(nms)) nms <- rep("", ...length())
  ifelse(nzchar(nms), nms, "<unnamed>")
}

rw_subset_outcome <- function(x, outcome) {
  if (is.null(outcome)) {
    return(x)
  }
  if (!is.character(outcome) || !length(outcome) || anyNA(outcome)) {
    rw_stop("`outcome` must be a character vector of outcome names.")
  }
  unknown <- setdiff(outcome, x$r2$outcome)
  if (length(unknown)) {
    rw_stop("Outcome(s) ", rw_collapse(unknown), " are not part of the ",
            "analysis. Available outcomes: ",
            rw_collapse(unique(x$r2$outcome)), ".")
  }
  x$weights <- x$weights[x$weights$outcome %in% outcome, , drop = FALSE]
  x$r2 <- x$r2[x$r2$outcome %in% outcome, , drop = FALSE]
  x
}

rw_has_intervals <- function(x, ci) {
  ci && !is.null(rw_interval_info(x))
}

# A factor identifying the regression of each row of x$weights, labelled
# with the regression and (optionally) its R-squared.
rw_regression_factor <- function(x, with_r2 = TRUE) {
  r2 <- x$r2
  lab <- rw_label(r2$outcome, r2$group, r2$level)
  if (with_r2) {
    lab <- paste0(lab, " (R\u00b2 = ", rw_fmt(r2$r2, 2), ")")
  }
  key_r2 <- paste(r2$block, r2$outcome, sep = "\r")
  key_w <- paste(x$weights$block, x$weights$outcome, sep = "\r")
  factor(lab[match(key_w, key_r2)], levels = lab)
}

# Order items by value within panels (largest value at the top of the axis).
rw_within <- function(item, value, panel) {
  key <- paste(item, as.integer(factor(panel)), sep = "\r")
  stats::reorder(factor(key), value)
}

rw_strip_within <- function(labels) sub("\r.*$", "", labels)

# Weights on the requested scale, with signs and interval bounds.
rw_weight_values <- function(x, scale, sign) {
  d <- x$weights
  value <- d[[scale]]
  if (is.null(rw_interval_info(x))) {
    lower <- upper <- rep(NA_real_, nrow(d))
  } else {
    lower <- d[[paste0(scale, "_lower")]]
    upper <- d[[paste0(scale, "_upper")]]
  }
  s <- switch(sign,
    none = rep(1L, nrow(d)),
    beta = d$sign_beta,
    orthogonal = d$sign_orth
  )
  s[s == 0L] <- 1L
  d$value <- value * s
  d$lower <- ifelse(s < 0L, -upper, lower)
  d$upper <- ifelse(s < 0L, -lower, upper)
  d$direction <- factor(ifelse(s < 0L, "Negative", "Positive"),
                        levels = c("Positive", "Negative"))
  d
}

rw_value_title <- function(scale, sign) {
  title <- if (scale == "rescaled") {
    "rescaled relative weight (% of R\u00b2)"
  } else {
    "raw relative weight"
  }
  if (sign == "none") {
    paste0(toupper(substring(title, 1, 1)), substring(title, 2))
  } else {
    paste("Signed", title)
  }
}

rw_caption <- function(x, sign, show_ci, what = "Error bars") {
  parts <- character()
  if (!is.null(x$draws)) {
    parts <- c(parts, "Estimates: posterior means.")
  }
  if (sign == "beta") {
    parts <- c(parts, "Signs: standardized coefficients.")
  } else if (sign == "orthogonal") {
    parts <- c(parts, "Signs: orthogonal-counterpart coefficients.")
  }
  if (show_ci) {
    iv <- rw_interval_info(x)
    unit <- if (iv$type == "credible") " posterior draws" else " replicates"
    parts <- c(parts, paste0(what, ": ", iv$text, " (", iv$n, unit, ")."))
  }
  if (length(parts)) paste(parts, collapse = "\n") else NULL
}

rw_scale_discrete <- function(aesthetic, n, name, values_order = seq_len(n)) {
  fun <- if (aesthetic == "fill") {
    ggplot2::scale_fill_manual
  } else {
    ggplot2::scale_colour_manual
  }
  if (n <= length(rw_palette)) {
    fun(values = rw_palette[values_order], name = name)
  } else if (aesthetic == "fill") {
    ggplot2::scale_fill_discrete(name = name)
  } else {
    ggplot2::scale_colour_discrete(name = name)
  }
}

rw_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title.position = "plot",
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold", hjust = 0)
    )
}

# Plots -------------------------------------------------------------------

rw_plot_weights <- function(x, scale, sign, ci) {
  d <- rw_weight_values(x, scale, sign)
  show_ci <- rw_has_intervals(x, ci)
  signed <- sign != "none"
  d$panel <- rw_regression_factor(x)
  d$y <- rw_within(d$predictor, d$value, d$panel)
  d$label <- if (scale == "rescaled") {
    paste0(rw_fmt(d$value, 1L), "%")
  } else {
    rw_fmt(d$value, 3L)
  }
  d$label_x <- d$value
  if (show_ci) {
    d$label_x <- ifelse(d$value >= 0,
                        pmax(d$value, d$upper, na.rm = TRUE),
                        pmin(d$value, d$lower, na.rm = TRUE))
  }
  d$hjust <- ifelse(d$value >= 0, -0.15, 1.15)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$value, y = .data$y))
  if (signed) {
    p <- p +
      ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
      ggplot2::geom_col(ggplot2::aes(fill = .data$direction), width = 0.7) +
      ggplot2::scale_fill_manual(
        values = c(Positive = "#0072B2", Negative = "#D55E00"),
        name = "Sign", drop = TRUE
      )
  } else {
    p <- p + ggplot2::geom_col(fill = "#0072B2", width = 0.7)
  }
  if (show_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
      width = 0.25, colour = "grey20", orientation = "y", na.rm = TRUE
    )
  }
  p <- p +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$label_x, label = .data$label,
                   hjust = .data$hjust),
      size = 3.2, na.rm = TRUE
    ) +
    ggplot2::scale_y_discrete(labels = rw_strip_within) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(
      mult = c(if (any(d$value < 0, na.rm = TRUE)) 0.2 else 0, 0.2)
    ))

  subtitle <- NULL
  if (nlevels(d$panel) > 1L) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$panel), scales = "free_y")
  } else {
    subtitle <- paste0("Outcome: ", levels(d$panel))
  }
  p +
    ggplot2::labs(title = "Relative weights", subtitle = subtitle,
                  x = rw_value_title(scale, sign), y = NULL,
                  caption = rw_caption(x, sign, show_ci)) +
    rw_theme()
}

rw_plot_r2 <- function(x, ci) {
  d <- x$weights
  r2 <- x$r2
  show_ci <- rw_has_intervals(x, ci)
  lab <- rw_label(r2$outcome, r2$group, r2$level)
  key <- paste(d$block, d$outcome, sep = "\r")
  d$regression <- factor(lab[match(key, paste(r2$block, r2$outcome,
                                               sep = "\r"))],
                         levels = rev(lab))
  d$predictor <- factor(d$predictor, levels = unique(d$predictor))
  d$end <- stats::ave(d$raw, key, FUN = cumsum)
  d$mid <- d$end - d$raw / 2

  totals <- data.frame(
    regression = factor(lab, levels = rev(lab)),
    r2 = r2$r2,
    lower = if (show_ci) r2$r2_lower else NA_real_,
    upper = if (show_ci) r2$r2_upper else NA_real_
  )
  totals$label <- paste0("R\u00b2 = ", rw_fmt(totals$r2, 3L))

  # Label only the segments that are wide enough to hold their label.
  axis_max <- max(c(totals$r2, totals$upper), na.rm = TRUE)
  d$segment_label <- ifelse(d$raw >= 0.07 * axis_max,
                            paste0(rw_fmt(d$rescaled, 0L), "%"), "")

  p <- ggplot2::ggplot(d, ggplot2::aes(y = .data$regression)) +
    ggplot2::geom_tile(
      ggplot2::aes(x = .data$mid, width = .data$raw, fill = .data$predictor),
      height = 0.55, colour = "white", linewidth = 0.3
    ) +
    ggplot2::geom_text(ggplot2::aes(x = .data$mid,
                                    label = .data$segment_label),
                       size = 3)
  if (show_ci) {
    # Drawn just below each bar so that it does not cross the labels.
    p <- p + ggplot2::geom_errorbar(
      data = totals,
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper,
                   y = .data$regression),
      inherit.aes = FALSE, width = 0.12, orientation = "y", na.rm = TRUE,
      position = ggplot2::position_nudge(y = -0.36)
    )
  }
  p +
    ggplot2::geom_text(
      data = totals,
      ggplot2::aes(x = .data$r2, y = .data$regression, label = .data$label),
      inherit.aes = FALSE, hjust = -0.1, size = 3.2
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.25))
    ) +
    rw_scale_discrete("fill", nlevels(d$predictor), "Predictor") +
    ggplot2::labs(title = "Decomposition of R\u00b2",
                  subtitle = "Segments are raw relative weights; labels give their share of R\u00b2",
                  x = "R\u00b2 (sum of the raw weights)", y = NULL,
                  caption = rw_caption(x, "none", show_ci,
                                       what = "Error bars for R\u00b2")) +
    rw_theme()
}

rw_plot_compare <- function(x) {
  d <- x$weights
  d$panel <- rw_regression_factor(x)
  key <- paste(d$block, d$outcome, sep = "\r")
  share <- function(v) 100 * v / sum(v)
  measures <- c("Relative weight", "Squared standardized coefficient",
                "Squared correlation")
  n <- nrow(d)
  long <- data.frame(
    panel = rep(d$panel, 3L),
    y = rep(rw_within(d$predictor, d$rescaled, d$panel), 3L),
    measure = factor(rep(measures, each = n), levels = rev(measures)),
    value = c(d$rescaled, stats::ave(d$beta^2, key, FUN = share),
              stats::ave(d$r^2, key, FUN = share))
  )

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value, y = .data$y,
                                          fill = .data$measure)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                      width = 0.75) +
    ggplot2::scale_y_discrete(labels = rw_strip_within) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(
      mult = c(0, 0.05)
    )) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(rw_palette[c(1L, 2L, 3L)], measures),
      breaks = measures, name = NULL
    )
  subtitle <- "Each index as a percentage of its total over the predictors"
  if (nlevels(d$panel) > 1L) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$panel), scales = "free_y")
  } else {
    subtitle <- paste0(subtitle, "; outcome: ", levels(d$panel))
  }
  p +
    ggplot2::labs(title = "Importance indices", subtitle = subtitle,
                  x = "Share of the total (%)", y = NULL) +
    rw_theme()
}

rw_plot_groups <- function(x, scale, sign, ci) {
  ngroups <- x$info$ngroups
  nlevels <- x$info$nlevels
  if (ngroups * nlevels < 2L) {
    rw_stop("type = \"groups\" compares weights across the groups or levels ",
            "of a lavaan or blavaan model; this analysis has a single group ",
            "and level.")
  }
  d <- rw_weight_values(x, scale, sign)
  show_ci <- rw_has_intervals(x, ci)
  signed <- sign != "none"
  if (ngroups > 1L && nlevels > 1L) {
    d$block_label <- paste(d$group, d$level, sep = ", ")
    legend <- "Group, level"
    title <- "Relative weights by group and level"
  } else if (ngroups > 1L) {
    d$block_label <- d$group
    legend <- "Group"
    title <- "Relative weights by group"
  } else {
    d$block_label <- d$level
    legend <- "Level"
    title <- "Relative weights by level"
  }
  # Reversed levels put the first group on top of each dodged set; the
  # legend is reversed back below.
  d$block_label <- factor(d$block_label, levels = rev(unique(d$block_label)))
  d$panel <- factor(d$outcome, levels = unique(d$outcome))
  d$y <- rw_within(d$predictor, d$value, d$panel)
  if (!show_ci) {
    d$lower <- d$upper <- d$value
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$value, xmin = .data$lower,
                                       xmax = .data$upper, y = .data$y,
                                       colour = .data$block_label))
  if (signed) {
    p <- p + ggplot2::geom_vline(xintercept = 0, colour = "grey40")
  }
  p <- p +
    ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.6),
                             orientation = "y", size = 0.4, na.rm = TRUE) +
    ggplot2::expand_limits(x = 0) +
    ggplot2::scale_y_discrete(labels = rw_strip_within) +
    rw_scale_discrete("colour", nlevels(d$block_label), legend,
                      values_order = rev(seq_len(nlevels(d$block_label)))) +
    ggplot2::guides(colour = ggplot2::guide_legend(reverse = TRUE))
  subtitle <- NULL
  if (nlevels(d$panel) > 1L) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$panel), scales = "free_y")
  } else {
    subtitle <- paste0("Outcome: ", levels(d$panel))
  }
  p +
    ggplot2::labs(title = title,
                  subtitle = subtitle, x = rw_value_title(scale, sign),
                  y = NULL, caption = rw_caption(x, sign, show_ci)) +
    rw_theme()
}
