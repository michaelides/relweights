#' relweights: relative weight analysis for linear models and structural
#' equation models
#'
#' Relative weight analysis (Johnson, 2000) divides the variance a regression
#' explains (\eqn{R^2}) among correlated predictors. Each weight combines a
#' predictor's unique contribution with the contribution it shares with the
#' other predictors, and the weights sum to \eqn{R^2}.
#'
#' The package works directly on fitted models:
#'
#' * [relweights()] computes the weights for a linear model fitted with
#'   [stats::lm()], or for every structural regression in a model fitted with
#'   'lavaan' (see [lavaan::sem()]). Optional bootstrap percentile intervals
#'   cover the raw weights, the rescaled weights and \eqn{R^2}.
#' * For Bayesian models fitted with 'brms' (gaussian linear models) or
#'   'blavaan' (structural equation models), [relweights()] computes the
#'   weights in every posterior draw and reports posterior means, credible
#'   intervals and posterior probabilities of signs and ranks.
#' * [autoplot.relweights()] and [plot.relweights()] draw the results with
#'   'ggplot2'.
#' * [print.relweights()] and [as.data.frame.relweights()] display and
#'   extract them.
#'
#' Start with `vignette("relweights", package = "relweights")` for a worked
#' introduction, and see `vignette("bayesian", package = "relweights")` for
#' Bayesian models.
#'
#' @references
#' Johnson, J. W. (2000). A heuristic method for estimating the relative weight
#' of predictor variables in multiple regression. *Multivariate Behavioral
#' Research, 35*(1), 1--19. \doi{10.1207/S15327906MBR3501_1}
#'
#' Johnson, J. W., & LeBreton, J. M. (2004). History and use of relative
#' importance indices in organizational research. *Organizational Research
#' Methods, 7*(3), 238--257. \doi{10.1177/1094428104266510}
#'
#' Tonidandel, S., & LeBreton, J. M. (2011). Relative importance analysis: A
#' useful supplement to regression analysis. *Journal of Business and
#' Psychology, 26*(1), 1--9. \doi{10.1007/s10869-010-9204-3}
#'
#' Tonidandel, S., & LeBreton, J. M. (2015). RWA Web: A free, comprehensive,
#' web-based, and user-friendly tool for relative weight analyses. *Journal of
#' Business and Psychology, 30*(2), 207--216. \doi{10.1007/s10869-014-9351-z}
#'
#' Tonidandel, S., LeBreton, J. M., & Johnson, J. W. (2009). Determining the
#' statistical significance of relative weights. *Psychological Methods,
#' 14*(4), 387--399. \doi{10.1037/a0017735}
#'
#' Thomas, D. R., Zumbo, B. D., Kwan, E., & Schweitzer, L. (2014). On
#' Johnson's (2000) relative weights method for assessing variable importance:
#' A reanalysis. *Multivariate Behavioral Research, 49*(4), 329--338.
#' \doi{10.1080/00273171.2014.905766}
#'
#' Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
#' *Journal of Statistical Software, 48*(2), 1--36. \doi{10.18637/jss.v048.i02}
#'
#' Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2019). R-squared for
#' Bayesian regression models. *The American Statistician, 73*(3), 307--309.
#' \doi{10.1080/00031305.2018.1549100}
#'
#' Bürkner, P.-C. (2017). brms: An R package for Bayesian multilevel models
#' using Stan. *Journal of Statistical Software, 80*(1), 1--28.
#' \doi{10.18637/jss.v080.i01}
#'
#' Merkle, E. C., Fitzsimmons, E., Uanhoro, J., & Goodrich, B. (2021).
#' Efficient Bayesian structural equation modeling in Stan. *Journal of
#' Statistical Software, 100*(6), 1--22. \doi{10.18637/jss.v100.i06}
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom ggplot2 autoplot .data
NULL
