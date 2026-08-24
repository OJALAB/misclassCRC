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
#' @param control An optional control object reserved for future use.
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
  control = NULL
) {

  call <- match.call()

  outcome_dist <- match.arg(
    outcome_dist,
    c("ztnegbin", "ztpois")
  )

  validate_crc_input(
    data = data,
    captures = captures,
    capture_formula = capture_formula,
    outcome = outcome,
    outcome_formula = outcome_formula,
    outcome_dist = outcome_dist,
    misclass = misclass,
    latent_classes = latent_classes,
    control = control
  )

  parsed_captures <- parse_captures(
    captures = captures,
    data = data
  )

  parsed_outcome <- parse_outcome(
    outcome = outcome,
    data = data
  )

  parsed_misclass <- parse_misclass(
    misclass = misclass,
    data = data,
    capture_names = parsed_captures$names
  )

  parsed_capture_formula <- parse_capture_formula(
    capture_formula = capture_formula,
    data = data,
    capture_names = parsed_captures$names,
    latent_classes = latent_classes
  )

  parsed_outcome_formula <- parse_outcome_formula(
    outcome_formula = outcome_formula,
    data = data,
    outcome = parsed_outcome,
    latent_classes = latent_classes
  )

  validate_model_compatibility(
    capture_model = parsed_capture_formula,
    outcome_model = parsed_outcome_formula
  )

  validate_misclass_compatibility(
    misclass = parsed_misclass,
    capture_model = parsed_capture_formula,
    outcome_model = parsed_outcome_formula
  )

  model_matrices <- build_model_matrices(
    data = data,
    captures = parsed_captures,
    capture_model = parsed_capture_formula,
    outcome = parsed_outcome,
    outcome_model = parsed_outcome_formula,
    misclass = parsed_misclass,
    latent_classes = latent_classes
  )

  structure(
    list(
      call = call,
      data = data,
      model_matrices = model_matrices,
      captures = parsed_captures,
      capture_model = parsed_capture_formula,
      outcome = parsed_outcome,
      outcome_model = parsed_outcome_formula,
      outcome_dist = outcome_dist,
      misclass = parsed_misclass,
      latent_classes = latent_classes,
      control = control,
      fitted = FALSE
    ),
    class = c("crcfit_unfitted", "crcfit")
  )

}