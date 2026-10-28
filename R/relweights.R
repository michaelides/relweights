#' Relative weight analysis
#'
#' `relweights()` computes Johnson's (2000) relative weights for a linear model
#' fitted with [stats::lm()] or for the structural regressions of a model
#' fitted with 'lavaan'. Relative weights divide the explained variance
#' (\eqn{R^2}) of a regression among its predictors, including the variance
#' that correlated predictors explain jointly.
#'
#' @section The method:
#' Let \eqn{R_{xx}} be the correlation matrix of the \eqn{p} predictors and
#' \eqn{r_{xy}} their correlations with the outcome. With the
#' eigendecomposition \eqn{R_{xx} = V D V'}, the matrix
#' \eqn{\Lambda = V D^{1/2} V'} holds the correlations between the predictors
#' and their *orthogonal counterparts* \eqn{Z}, the uncorrelated variables
#' that are closest to the predictors in the least-squares sense. Regressing
#' the outcome on \eqn{Z} gives the coefficients
#' \eqn{\gamma = \Lambda^{-1} r_{xy}} (column `orth`), and because the columns
#' of \eqn{Z} are uncorrelated, \eqn{R^2 = \sum_k \gamma_k^2}. The variance
#' explained by each \eqn{Z_k} is shared out among the predictors in
#' proportion to \eqn{\lambda_{jk}^2}, giving the raw weight
#' \deqn{\epsilon_j = \sum_k \lambda_{jk}^2 \gamma_k^2.}
#' Raw weights are non-negative and sum to \eqn{R^2}; rescaled weights express
#' them as percentages of \eqn{R^2} and sum to 100.
#'
#' Relative weights approximate the general dominance weights of dominance
#' analysis at a fraction of the computing cost, but they are a heuristic:
#' see Thomas et al. (2014) for a critique of their derivation.
#'
#' @section Signs:
#' Relative weights carry no sign. Two signs are reported so that the
#' direction of each association can be shown alongside its importance:
#' \describe{
#'   \item{`sign_beta`}{the sign of the standardized coefficient of the
#'     predictor in the full regression (`beta`), which agrees with the
#'     model's own coefficient;}
#'   \item{`sign_orth`}{the sign of the coefficient of the predictor's
#'     orthogonal counterpart (`orth`), the convention of RWA Web and of the
#'     'rwa' package.}
#' }
#' The two usually agree but can differ when predictors are strongly
#' correlated (for example under suppression). The zero-order correlation is
#' reported in column `r`.
#'
#' @section Linear models:
#' The predictors are the columns of the model matrix without the intercept,
#' named as in `coef(object)`, and the outcome is the response. Correlations
#' are computed from the rows used in the fit. With prior weights, weighted
#' correlations are used (rows with zero weight are dropped), which reproduces
#' the weighted \eqn{R^2} of the model.
#'
#' A term that produces several model-matrix columns, such as a factor with
#' more than two levels, [stats::poly()], a spline basis or a matrix, receives one
#' weight per column and a warning. For factors, polynomials and splines these
#' weights depend on the chosen contrasts or basis, and so do the weights of
#' the other predictors: changing the reference level of a factor changes all
#' weights, although \eqn{R^2} does not change. Interaction and power columns
#' are correlated with their components unless the components are centred,
#' which also affects the weights.
#'
#' The model must have an intercept, no offset and no aliased coefficients,
#' and its class must be exactly `"lm"`.
#'
#' @section lavaan models:
#' Every regression defined with `~` is analysed, separately in every block
#' (group and level) of the model. The correlations are the model-implied
#' correlations among all observed and latent variables
#' (`lavaan::lavInspect(object, "cor.all")`), so regressions on latent
#' variables are analysed free of the measurement error captured by the
#' measurement model. When the residual of an outcome is uncorrelated with its
#' predictors, as in most models, the regression reproduces the standardized
#' coefficients and the \eqn{R^2} reported by 'lavaan'. Otherwise (for
#' example with a residual covariance between the outcome and one of its
#' predictors, or in non-recursive models) the weights decompose the
#' model-implied linear regression of the outcome on its predictors instead,
#' and a warning is given. For ordered categorical variables, correlations
#' refer to the underlying latent responses.
#'
#' Exogenous covariates are excluded from the implied correlation matrix when
#' 'lavaan' conditions on them (`conditional.x = TRUE`, the default when some
#' endogenous variables are ordered); refit with `conditional.x = FALSE` to
#' analyse them. Column `n` of the `r2` table counts observations at the first
#' level and clusters at the second. The model must have converged, and a
#' warning is given when its solution is not admissible.
#'
#' @section Bootstrap intervals:
#' For linear models, rows of the data (with their prior weights) are
#' resampled with [boot::boot()] and the analysis is repeated in each
#' resample. For 'lavaan' models, the model is refitted to each resample with
#' [lavaan::lavBootstrap()] and the implied correlations are recomputed; this
#' needs a model fitted to raw data and is not available for multilevel
#' models. The intervals are percentile intervals; replicates that cannot be
#' analysed are excluded with a warning. Use [set.seed()] beforehand for
#' reproducible intervals, and at least 1000 replicates for reporting.
#'
#' Raw weights are non-negative, so an interval that excludes zero does not
#' show that a weight differs from zero (Tonidandel et al., 2009). Rescaled
#' weights are proportions of \eqn{R^2} that sum to 100, so their intervals
#' are best read descriptively.
#'
#' @param object A fitted model: a linear model from [stats::lm()] or a model
#'   fitted with 'lavaan', for example with [lavaan::sem()].
#' @param outcome For 'lavaan' models, an optional character vector naming the
#'   outcomes (left-hand sides of `~` regressions) to analyse. By default every
#'   regression is analysed.
#' @param bootstrap Logical. Compute bootstrap percentile confidence intervals
#'   for the raw weights, the rescaled weights and \eqn{R^2}?
#' @param n_boot Number of bootstrap replicates.
#' @param conf_level Confidence level of the intervals.
#' @param ... Options for the bootstrap: `parallel`, `ncpus` and `cl` are
#'   passed to [boot::boot()] for linear models and to
#'   [lavaan::lavBootstrap()] for 'lavaan' models; `iseed` and `type` are also
#'   accepted for 'lavaan' models. Any other argument is an error.
#'
#' @return An object of class `relweights`: a list with components
#' \describe{
#'   \item{`weights`}{A data frame with one row per predictor of each
#'     regression: `block`, `group` and `level` (the 'lavaan' block, group
#'     label and level label, `NA` when not applicable), `outcome`,
#'     `predictor`, `term` (the model term a column belongs to), `raw` (raw
#'     weights, summing to \eqn{R^2}), `rescaled` (percentages of
#'     \eqn{R^2}), `beta` (standardized coefficients), `orth` (coefficients
#'     of the orthogonal counterparts), `r` (zero-order correlations with the
#'     outcome), and `sign_beta` and `sign_orth` (the signs of `beta` and
#'     `orth`, coded -1, 0 and 1). With `bootstrap = TRUE`, also `raw_lower`,
#'     `raw_upper`, `rescaled_lower` and `rescaled_upper`.}
#'   \item{`r2`}{A data frame with one row per regression: `block`, `group`,
#'     `level`, `outcome`, `n`, `r2` (the sum of the raw weights) and
#'     `r2_model` (the \eqn{R^2} reported by the fitted model), plus
#'     `r2_lower` and `r2_upper` with `bootstrap = TRUE`.}
#'   \item{`matrices`}{A list with one element per regression holding `R`,
#'     the correlation matrix of the outcome (first) and the predictors, and
#'     `lambda`, the matrix \eqn{\Lambda}.}
#'   \item{`boot`}{`NULL` without bootstrap; otherwise a list with the
#'     replicate matrix `t` (one row per replicate, `NA` rows for failed
#'     replicates), the estimates `t0`, an `index` data frame describing the
#'     columns of `t`, `n_boot`, `n_failed`, `n_nonadmissible`, `conf_level`
#'     and `method`.}
#'   \item{`info`}{A list describing the model: `model_class`, `call`,
#'     `weighted` (whether a linear model had prior weights; `NA` for
#'     'lavaan' models), `ngroups`, `nlevels`, `group_labels`,
#'     `level_labels` and `multi_column_terms` (the terms of a linear model
#'     that produced several model-matrix columns).}
#' }
#'
#' @references
#' Johnson, J. W. (2000). A heuristic method for estimating the relative weight
#' of predictor variables in multiple regression. *Multivariate Behavioral
#' Research, 35*(1), 1--19. \doi{10.1207/S15327906MBR3501_1}
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
#' @seealso [autoplot.relweights()] for plots, [print.relweights()] and
#'   [as.data.frame.relweights()].
#'
#' @examples
#' # A linear model
#' fit <- lm(mpg ~ cyl + disp + hp + gear, data = mtcars)
#' rw <- relweights(fit)
#' rw
#' as.data.frame(rw)
#'
#' # A regression among latent variables
#' model <- "
#'   ind60 =~ x1 + x2 + x3
#'   dem60 =~ y1 + y2 + y3 + y4
#'   dem65 =~ y5 + y6 + y7 + y8
#'   dem65 ~ ind60 + dem60
#' "
#' sem_fit <- lavaan::sem(model, data = lavaan::PoliticalDemocracy)
#' relweights(sem_fit)
#'
#' # Latent and observed predictors, analysed in two groups
#' hs_model <- "
#'   visual  =~ x1 + x2 + x3
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#'   speed ~ visual + textual + ageyr
#' "
#' hs_fit <- lavaan::sem(hs_model, data = lavaan::HolzingerSwineford1939,
#'                       group = "school")
#' relweights(hs_fit)
#'
#' \donttest{
#' # Bootstrap percentile intervals
#' set.seed(2024)
#' relweights(fit, bootstrap = TRUE, n_boot = 500)
#' }
#' @export
relweights <- function(object, ...) {
  UseMethod("relweights")
}

#' @rdname relweights
#' @export
relweights.default <- function(object, ...) {
  rw_stop("relweights() supports linear models fitted with stats::lm() and ",
          "models fitted with lavaan; objects of class ",
          rw_collapse(class(object)[1L]), " are not supported.")
}
