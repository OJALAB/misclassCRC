#' Prepare Capture-Recapture Model Inputs
#' 
#' @noRd
prepare_crc <- function(
  data,
  captures,
  capture_formula,
  outcome,
  outcome_formula,
  outcome_dist,
  misclass,
  latent_classes,
  control,
  verbose
) {

  outcome_dist <- match.arg(
    outcome_dist,
    c("ztnegbin", "ztpois")
  )

  control <- parse_control(control)

  validate_crc_input(
    data = data,
    captures = captures,
    capture_formula = capture_formula,
    outcome = outcome,
    outcome_formula = outcome_formula,
    outcome_dist = outcome_dist,
    misclass = misclass,
    latent_classes = latent_classes,
    verbose = verbose
  )

  captures <- parse_captures(
    captures = captures,
    data = data
  )

  outcome <- parse_outcome(
    outcome = outcome,
    data = data
  )

  misclass <- parse_misclass(
    misclass = misclass,
    data = data,
    capture_names = captures$names
  )

  capture_model <- parse_capture_formula(
    capture_formula = capture_formula,
    data = data,
    capture_names = captures$names,
    latent_classes = latent_classes
  )

  outcome_model <- parse_outcome_formula(
    outcome_formula = outcome_formula,
    data = data,
    outcome = outcome,
    latent_classes = latent_classes
  )

  validate_model_compatibility(
    capture_model = capture_model,
    outcome_model = outcome_model
  )

  validate_misclass_compatibility(
    misclass = misclass,
    capture_model = capture_model,
    outcome_model = outcome_model
  )

  model_matrices <- build_model_matrices(
    data = data,
    captures = captures,
    capture_model = capture_model,
    outcome = outcome,
    outcome_model = outcome_model,
    misclass = misclass,
    latent_classes = latent_classes
  )

  initialization <- initialize_crc(
    model_matrices = model_matrices,
    misclass = misclass,
    latent_classes = latent_classes,
    control = control
  )

  structure(
    list(
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
      verbose = verbose
    ),
    class = c("crcfit_preparation")
  )

}