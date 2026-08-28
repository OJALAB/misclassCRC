#' Fit a Capture-Recapture Model
#' 
#' @param data A `data.frame` or `data.table` with one row per observed unit.
#' @param captures A one-sided formula identifying binary capture indicators, e.g.,
#' `~ source_1 + source_2 + source_3`.
#' @param capture_formula A one-sided formula specifying the log-linear capture model.
#' The special term `.latent` denotes the latent class.
#' @param outcome Either `NULL`, a single column name for an exactly observed outcome,
#' or a named character vector `c(lower = "lower_column", upper = "upper_column")`.
#' @param outcome_formula A one-sided formula for the mean of the outcome distribution.
#' Must be `NULL` when `outcome = NULL`.
#' @param outcome_dist The zero-truncated outcome distribution.
#' Currently, `"ztnegbin"` and `"ztpois"` are supported.
#' @param misclass A specification of misclassification mechanisms.
#' @param latent_classes A positive integer giving the number of latent classes.
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
#' be displayed during model fitting. The default is `FALSE`.
#'
#' @return 
#' An object of class `"crcfit"`.
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

  structure(
    list(
      call = call,
      data = data,
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
      fitted = FALSE
    ),
    class = c("crcfit_unfitted", "crcfit")
  )

}