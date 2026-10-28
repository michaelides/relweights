#' Relative weight analysis
#'
#' `relweights()` computes Johnson's (2000) relative weights for a linear model
#' fitted with [stats::lm()], for the structural regressions of a model
#' fitted with 'lavaan' or 'blavaan', or for a linear model fitted with
#' 'brms'. Relative weights divide the explained variance (\eqn{R^2}) of a
#' regression among its predictors, including the variance that correlated
#' predictors explain jointly. For the Bayesian models fitted with 'blavaan'
#' and 'brms', the weights are computed in every posterior draw.
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
#' @section Bayesian models:
#' For models fitted with 'blavaan' or 'brms', the analysis is repeated in
#' every posterior draw, which gives the posterior distribution of the
#' weights. The main columns hold posterior means, which keep the weights
#' additive: the raw weights sum to the posterior mean of \eqn{R^2} and the
#' rescaled weights to 100. The interval columns hold equal-tailed credible
#' intervals. Further columns give the posterior probabilities that `beta`
#' and `orth` are positive, and the R-hat and bulk effective sample size of
#' each raw weight (computed with the 'posterior' package); a warning is
#' given when an R-hat exceeds 1.05. With `what = "ranks"` or
#' `what = "pairwise"`, [as.data.frame.relweights()] returns the posterior
#' probabilities that each predictor has a given rank, or a larger weight
#' than each other predictor. Like bootstrap intervals, credible intervals of
#' weights never include zero, so they do not test whether a weight differs
#' from zero.
#'
#' **'brms' models.** Supported models use the gaussian family with the
#' identity link and population-level effects only. Multivariate models
#' (`mvbind()`) are analysed one response at a time. Group-level terms,
#' predicted distributional parameters (such as `sigma ~ x`), non-linear
#' formulas, smooth, special, Gaussian process and autocorrelation terms,
#' offsets and addition terms (`y | ...`) are errors. The predictors are the
#' columns of the design matrix without the intercept (see
#' [brms::standata()]), and terms with several columns are treated as for
#' linear models. Because 'brms' conditions on the predictors, their
#' correlations are the sample correlations in every draw; each draw
#' contributes its regression coefficients \eqn{b} and residual standard
#' deviation \eqn{\sigma}. With the covariance matrix \eqn{S} of the
#' predictors, a draw implies the outcome variance \eqn{b'Sb + \sigma^2}, the
#' standardized coefficients and correlations, and the model-based
#' \eqn{R^2 = b'Sb / (b'Sb + \sigma^2)} of Gelman et al. (2019), which the
#' raw weights sum to. Column `r2_model` holds the posterior mean of the
#' residual-based \eqn{R^2} of [brms::bayes_R2()], which uses the variance of
#' the residuals in place of \eqn{\sigma^2}. The intervals reflect the
#' uncertainty of \eqn{b} and \eqn{\sigma}, not that of the predictor
#' correlations, which the bootstrap of linear models also covers.
#'
#' **'blavaan' models.** In each draw, the parameters are put into the model,
#' the model-implied correlations among all observed and latent variables
#' are recomputed in every group and level, and every regression is analysed
#' as for 'lavaan' models; multilevel models are supported. The intervals
#' therefore reflect the uncertainty of all model parameters, except for the
#' covariances of exogenous observed covariates, which are fixed at their
#' sample values when `fixed.x = TRUE` (the default). Column `r2_model` holds
#' the posterior mean of the model's own \eqn{R^2}, computed from the
#' residual variance in each draw. Draws that cannot be analysed (for
#' example, because an implied correlation matrix is not positive
#' semi-definite) are excluded with a warning.
#'
#' @param object A fitted model: a linear model from [stats::lm()], a model
#'   fitted with 'lavaan' (for example with [lavaan::sem()]) or 'blavaan'
#'   (for example with [blavaan::bsem()]), or a model fitted with
#'   [brms::brm()].
#' @param outcome For 'lavaan', 'blavaan' and multivariate 'brms' models, an
#'   optional character vector naming the outcomes to analyse (left-hand
#'   sides of `~` regressions, or responses). By default every regression is
#'   analysed.
#' @param bootstrap Logical. Compute bootstrap percentile confidence intervals
#'   for the raw weights, the rescaled weights and \eqn{R^2}?
#' @param n_boot Number of bootstrap replicates.
#' @param conf_level Confidence level of the bootstrap intervals, or
#'   probability of the credible intervals for Bayesian models.
#' @param ndraws For Bayesian models, the maximum number of posterior draws to
#'   analyse. Evenly spaced draws are taken from every chain (at least one per
#'   chain). By default, all draws are analysed.
#' @param ... Options for the bootstrap: `parallel`, `ncpus` and `cl` are
#'   passed to [boot::boot()] for linear models and to
#'   [lavaan::lavBootstrap()] for 'lavaan' models; `iseed` and `type` are also
#'   accepted for 'lavaan' models. Any other argument is an error, and for
#'   Bayesian models `...` must be empty.
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
#'     `raw_upper`, `rescaled_lower` and `rescaled_upper`. For Bayesian
#'     models, the estimates are posterior means, the interval columns hold
#'     credible intervals, and the further columns `p_beta_pos` and
#'     `p_orth_pos` (posterior probabilities that `beta` and `orth` are
#'     positive) and `rhat` and `ess_bulk` (R-hat and bulk effective sample
#'     size of the raw weight) are added.}
#'   \item{`r2`}{A data frame with one row per regression: `block`, `group`,
#'     `level`, `outcome`, `n`, `r2` (the sum of the raw weights) and
#'     `r2_model` (the \eqn{R^2} reported by the fitted model; for Bayesian
#'     models, see the section on Bayesian models), plus `r2_lower` and
#'     `r2_upper` with bootstrap or credible intervals.}
#'   \item{`matrices`}{A list with one element per regression holding `R`,
#'     the correlation matrix of the outcome (first) and the predictors, and
#'     `lambda`, the matrix \eqn{\Lambda}.}
#'   \item{`boot`}{`NULL` without bootstrap; otherwise a list with the
#'     replicate matrix `t` (one row per replicate, `NA` rows for failed
#'     replicates), the estimates `t0`, an `index` data frame describing the
#'     columns of `t`, `n_boot`, `n_failed`, `n_nonadmissible`, `conf_level`
#'     and `method`.}
#'   \item{`draws`}{`NULL` except for Bayesian models. Otherwise a list with
#'     the `chain` and `iteration` of every analysed draw, `ok` (whether the
#'     draw could be analysed), `conf_level`, `n_failed` and `regressions`:
#'     one element per regression, named as in `matrices`, holding the
#'     per-draw matrices `raw`, `rescaled`, `beta`, `orth` and `r` (draws by
#'     predictors) and the vectors `r2` and `r2_model`.}
#'   \item{`info`}{A list describing the model: `model_class`, `call`,
#'     `weighted` (whether a linear model had prior weights), `ngroups`,
#'     `nlevels`, `group_labels`, `level_labels` and `multi_column_terms`
#'     (the terms that produced several model-matrix columns). For Bayesian
#'     models, also `intervals` (`"credible"`), `ndraws`, `nchains` and
#'     `r2_model_label`, and for 'brms' models `formula`.}
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
#' Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2019). R-squared for
#' Bayesian regression models. *The American Statistician, 73*(3), 307--309.
#' \doi{10.1080/00031305.2018.1549100}
#'
#' @seealso [autoplot.relweights()] for plots, [print.relweights()] and
#'   [as.data.frame.relweights()]; `vignette("bayesian", package =
#'   "relweights")` for Bayesian models.
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
#'
#' # A Bayesian structural equation model
#' if (requireNamespace("blavaan", quietly = TRUE)) {
#'   library(blavaan)
#'   bsem_fit <- bsem("x9 ~ x7 + x8 + ageyr",
#'                    data = lavaan::HolzingerSwineford1939,
#'                    burnin = 500, sample = 500, n.chains = 2, seed = 1,
#'                    bcontrol = list(refresh = 0))
#'   relweights(bsem_fit)
#' }
#' }
#'
#' \dontrun{
#' # A Bayesian linear model (compiles a Stan model)
#' brms_fit <- brms::brm(mpg ~ cyl + disp + hp + gear, data = mtcars,
#'                       seed = 1)
#' relweights(brms_fit)
#' }
#' @export
relweights <- function(object, ...) {
  UseMethod("relweights")
}

#' @rdname relweights
#' @export
relweights.default <- function(object, ...) {
  rw_stop("relweights() supports linear models fitted with stats::lm(), ",
          "models fitted with lavaan or blavaan, and brms models; objects of ",
          "class ", rw_collapse(class(object)[1L]), " are not supported.")
}
