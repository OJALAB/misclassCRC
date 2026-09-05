#' Summarize a Capture-Recapture Model Fit
#'
#' @param object An object of class `"crc_fit"`.
#' @param ... Additional arguments. Currently unused.
#' @return An object of class `"summary.crc_fit"`.
#' @export
summary.crc_fit <- function(
  object,
  ...
) {

  misclassified <- if (is.null(object$misclass)) {
    character()
  } else {
    vapply(
      object$misclass$specifications,
      function(specification) {
        truth <- if (length(specification$true_sources)) {
          paste0(
            "; true for ",
            paste(specification$true_sources, collapse = ", ")
          )
        } else {
          ""
        }
        paste0(length(specification$levels), " levels", truth)
      },
      character(1L)
    )
  }

  structure(
    list(
      call = object$call,
      model = list(
        "Observed units" = nrow(object$data),
        "Capture sources" = paste(object$captures$names, collapse = ", "),
        "Latent classes" = object$latent_classes,
        "Capture cells" = nrow(object$model_matrices$capture$cells),
        "Observation states" = nrow(object$model_matrices$states$data),
        "Misclassified variables" = if (length(misclassified)) {
          paste(names(misclassified), collapse = ", ")
        } else {
          "none"
        },
        "Outcome distribution" = if (is.null(object$outcome_model)) {
          "none"
        } else {
          switch(
            object$outcome_dist,
            ztnegbin = "zero-truncated negative binomial",
            ztpois = "zero-truncated Poisson"
          )
        }
      ),
      misclassification = misclassified,
      capture = list(
        formula = object$capture_model$formula,
        coefficients = object$em_fit$capture_fit$coefficients
      ),
      outcome = if (is.null(object$outcome_model)) {
        NULL
      } else {
        list(
          formula = object$outcome_model$formula,
          coefficients = object$em_fit$outcome_fit$coefficients,
          sigma = object$em_fit$outcome_fit$sigma
        )
      },
      convergence = object$em_fit$convergence,
      converged = object$em_fit$converged,
      iterations = object$em_fit$iterations
    ),
    class = "summary.crc_fit"
  )

}