#' Fit a Capture-Recapture Model
#'
#' Fits a log-linear capture-recapture model with optional latent classes,
#' misclassification correction for categorical covariates, and a
#' zero-truncated count outcome with missing and censored values.
#' The population size and the outcome total are then estimated with [predict.crc_fit()].
#'
#' @param data A `data.frame` or `data.table` with one row per observed unit,
#' i.e., each unit must be captured by at least one source.
#' @param captures A one-sided formula with at least two capture indicators
#' joined by `+`, e.g., `~ source_1 + source_2 + source_3`. The indicators
#' must be numeric or logical, coded as 0/1, and contain no missing values.
#' With two sources, a warning is issued because the number of unobserved
#' units is identifiable only under restrictive assumptions.
#' @param capture_formula A one-sided formula specifying the log-linear model
#' for the expected counts in the cells of the complete contingency table
#' cross-classified by capture histories, categorical covariates, and latent
#' classes. The following rules apply:
#'
#' - The formula must include an intercept and every capture indicator as a
#' main effect.
#' - Other variables must be factor or character columns without missing
#' values. Unused factor levels are dropped, and character variables are
#' converted to factors with sorted levels
#' (for misclassified variables, the levels follow the confusion matrix).
#' - Terms must be untransformed variable names or their interactions
#' (written with `:`), and every component of an interaction must also appear
#' as a main effect.
#' - Capture indicators may interact with each other, and each of them
#' individually with `.latent`, but not with other covariates.
#' - The special term `.latent` denotes the latent class. It must be included
#' as a main effect when `latent_classes > 1` and cannot be used otherwise.
#' `.latent` is a reserved name and cannot be a column of `data`.
#' @param outcome Either `NULL` (no outcome model), a single column name, or a
#' named character vector `c(lower = "lower_column", upper = "upper_column")`.
#' The outcome is a count variable taking values in \eqn{\{1, 2, \ldots\}}{{1, 2, ...}}.
#'
#' - If a single column name is given, the column must be numeric, its
#' non-missing values must be positive integers, and `NA` denotes a missing
#' outcome.
#' - If lower and upper bounds are given, both columns must be numeric, and
#' for each unit either both bounds are `NA` or both are observed. Lower
#' bounds must be positive integers, upper bounds must be positive integers
#' or `Inf`, and lower bounds must not exceed upper bounds. The status of
#' each outcome is determined as follows:
#'   - both bounds `NA` -- missing outcome,
#'   - `lower == upper` -- exactly observed outcome,
#'   - `upper == Inf` -- right-censored outcome, i.e., the outcome is at least
#'   `lower` (e.g., `lower = 2` and `upper = Inf` represent "2 or more"),
#'   - otherwise -- interval-censored outcome, i.e., the outcome lies between
#'   `lower` and `upper` (inclusive).
#'
#' Missing and censored outcomes are assumed to be missing or coarsened at
#' random given the variables in `outcome_formula`. At least one outcome must
#' be observed or censored.
#' @param outcome_formula A one-sided formula with an intercept for the
#' logarithm of the mean of the untruncated outcome distribution. It may
#' contain only variables from `capture_formula` other than the capture
#' indicators, `.latent` (when `latent_classes > 1`), and their interactions
#' (whose components must also appear as main effects). It must be supplied
#' when `outcome` is not `NULL` and must be `NULL` otherwise.
#' @param outcome_dist The zero-truncated outcome distribution: `"ztnegbin"`
#' (zero-truncated negative binomial, the default) or `"ztpois"`
#' (zero-truncated Poisson). For `"ztnegbin"`, the untruncated distribution
#' has mean \eqn{\mu} and variance \eqn{\mu + \sigma \mu^2}{mu + sigma * mu^2}, where
#' \eqn{\sigma > 0} is a dispersion parameter.
#' @param misclass `NULL` (no misclassification) or a named list with one
#' element per categorical variable subject to misclassification. Each name
#' must be a column of `data` used in `capture_formula` or `outcome_formula`,
#' and that column must be a factor or a character vector without missing
#' values. Each element is a list with the following components:
#'
#' - `matrix` -- a square confusion matrix whose row and column names contain
#' the same category levels in the same order. Rows correspond to true
#' categories and columns to observed categories, i.e., the element `[g, h]`
#' is the probability of observing category `h` when the true category is
#' `g`. Entries must be finite and non-negative, and each row must sum to 1.
#' All observed levels of the variable must appear in the matrix. The row
#' names define the true categories used in the model.
#' - `true_if` (optional) -- a one-sided formula with capture variables joined
#' by `+`, e.g., `~ source_1`. Units captured by at least one of these
#' sources are treated as correctly classified, i.e., their observed category
#' is taken as the true category. The confusion matrix applies to all other
#' units. If `true_if` is omitted, it applies to all units.
#'
#' Confusion matrices are treated as known (e.g., estimated beforehand from an
#' independent validation sample), so their estimation uncertainty is not
#' propagated. For several variables, the misclassification mechanisms are
#' assumed to be conditionally independent given the true categories.
#' @param latent_classes A positive integer giving the number of latent
#' classes. The default `1L` corresponds to a model without latent classes.
#' For `latent_classes > 1`, initial class weights are random, so use
#' \link[base:Random]{set.seed()} for reproducible results.
#' @param control `NULL` or a named list controlling initialization and model fitting.
#' If `NULL`, the default values are used. Supported elements are:
#'
#' - `init_alpha`: A finite positive number controlling the concentration of the
#' symmetric Dirichlet distribution used to initialize latent-class weights.
#' The default is `20`.
#' - `em`: A named list with `max_iter`, the maximum number of EM iterations
#' (default `1000L`), and `tolerance`, the convergence tolerance (default `1e-6`).
#' - `capture`: A named list with `max_iter`, the maximum number of iterations
#' used in the capture-model M-step (default `100L`), and `tolerance`, its
#' convergence tolerance (default `1e-8`).
#' - `outcome`: A named list with `max_iter`, the maximum number of iterations
#' used in the outcome-model M-step (default `1000L`), and `relative_tolerance`,
#' its relative tolerance (default `1e-8`).
#'
#' Unknown control elements are rejected.
#' @param verbose A logical value indicating whether progress information should
#' be displayed during model fitting. The default is `FALSE`. Currently unsupported.
#'
#' @return
#' An object of class `"crc_fit"`, i.e., a list with the following components:
#'
#' - `call` -- the matched call,
#' - `data` -- the input data,
#' - `em_fit` -- the results of the EM algorithm: `capture_fit` and
#' `outcome_fit` (the final capture-model and outcome-model fits with their
#' `coefficients`; `outcome_fit` also contains the dispersion `sigma` for
#' `"ztnegbin"` and is `NULL` without an outcome model), `state_weights` (the
#' posterior state probabilities), `cell_counts` (the expected cell counts),
#' `iterations`, `converged`, and `convergence` (the final changes used in the
#' stopping rule),
#' - `model_matrices` -- the cell grid, design matrices, and expanded states,
#' - `initialization` -- the initial state weights and cell counts,
#' - `captures`, `capture_model`, `outcome`, `outcome_model`, and `misclass`
#' -- the parsed model specifications,
#' - `outcome_dist`, `latent_classes`, `control`, and `verbose` -- the fitting
#' settings, with `control` completed with default values,
#' - `fitted` -- `TRUE` for a fitted model.
#' 
#' @references
#' Baffour, B., Brown, J. J., and Smith, P. W. F. (2021). Latent class
#' analysis for estimating an unknown population size -- with application to
#' censuses. *Journal of Official Statistics*, *37*(3), 673--697.
#' \doi{10.2478/jos-2021-0030}
#'
#' \enc{Brännäs}{Brannas}, K. (1992). Limited dependent Poisson regression.
#' *Journal of the Royal Statistical Society. Series D (The Statistician)*,
#' *41*(4), 413--423. \doi{10.2307/2349006}
#'
#' Struzik, A. (2026). *Estimation of job vacancy statistics under
#' misclassification and missing data using capture-recapture*.
#' Royal Statistical Society 2026 International Conference, Bournemouth, United Kingdom.
#' \url{https://github.com/OJALAB/misclassCRC/blob/main/references/Struzik-RSS-2026.pdf}
#' 
#' @examples
#' set.seed(123)
#' N <- 1000
#' population <- data.frame(
#'   source_1 = rbinom(N, 1, 0.5),
#'   source_2 = rbinom(N, 1, 0.4),
#'   source_3 = rbinom(N, 1, 0.3),
#'   group = sample(c("a", "b"), N, replace = TRUE)
#' )
#' population$outcome <- 1 + rpois(N, ifelse(population$group == "a", 1, 3))
#'
#' # Keep only the units captured by at least one source
#' sources <- c("source_1", "source_2", "source_3")
#' data <- population[rowSums(population[, sources]) > 0, ]
#' data$outcome[sample(nrow(data), 50)] <- NA
#'
#' fit <- crc_fit(
#'   data,
#'   captures = ~ source_1 + source_2 + source_3,
#'   capture_formula = ~ source_1 + source_2 + source_3 + group,
#'   outcome = "outcome",
#'   outcome_formula = ~ group
#' )
#' fit
#' predict(fit, by = ~ group)
#'
#' @family fitting functions
#' 
#' @export
crc_fit <- function(
  data,
  captures,
  capture_formula,
  outcome = NULL,
  outcome_formula = NULL,
  outcome_dist = c("ztnegbin", "ztpois"),
  misclass = NULL,
  latent_classes = 1L,
  control = NULL,
  verbose = FALSE
) {

  call <- match.call()

  preparation <- prepare_crc(
    data = data,
    captures = captures,
    capture_formula = capture_formula,
    outcome = outcome,
    outcome_formula = outcome_formula,
    outcome_dist = outcome_dist,
    misclass = misclass,
    latent_classes = latent_classes,
    control = control,
    verbose = verbose
  )

  model_matrices <- NULL
  initialization <- NULL
  capture_model <- NULL
  outcome_model <- NULL

  list2env(preparation, envir = environment())
  rm(preparation)

  em_fit <- perform_em_crc(
    model_matrices = model_matrices,
    initialization = initialization,
    outcome_dist = outcome_dist,
    control = control,
    verbose = verbose
  )

  if (!em_fit$converged) {
    warning(
      "The EM algorithm did not converge after ",
      em_fit$iterations,
      " iterations.",
      call. = FALSE
    )
  }

  structure(
    list(
      call = call,
      data = data,
      em_fit = em_fit,
      model_matrices = model_matrices,
      initialization = initialization,
      captures = captures,
      capture_model = capture_model,
      outcome = outcome,
      outcome_model = outcome_model,
      outcome_dist = outcome_dist,
      misclass = misclass,
      latent_classes = latent_classes,
      control = control,
      verbose = verbose,
      fitted = TRUE
    ),
    class = "crc_fit"
  )

}
