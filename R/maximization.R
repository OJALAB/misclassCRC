#' Maximize the Capture-Model Objective
#' 
#' @importFrom stats glm.control glm.fit poisson
#' 
#' @noRd
maximize_capture_crc <- function(
  model_matrices,
  cell_counts,
  start = NULL,
  glm_control = list()
) {

  matrix <- model_matrices$capture$matrix
  n_parameters <- ncol(matrix)

  if (
    length(cell_counts) != nrow(matrix) ||
      anyNA(cell_counts) ||
      any(!is.finite(cell_counts)) ||
      any(cell_counts < 0)
  ) {
    stop("Expected capture-cell counts are invalid.", call. = FALSE)
  }

  if (is.null(start)) {
    start <- numeric(n_parameters)
  }

  if (
    length(start) != n_parameters ||
      anyNA(start) ||
      any(!is.finite(start))
  ) {
    stop("Capture-model starting values are invalid.", call. = FALSE)
  }

  fit <- withCallingHandlers(
    glm.fit(
      x = matrix,
      y = cell_counts,
      family = poisson(link = "log"),
      start = start,
      control = do.call(glm.control, glm_control)
    ),
    warning = function(warning) {
      if (grepl("non-integer x", conditionMessage(warning), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )

  coefficients <- fit$coefficients
  fitted_mean <- fit$fitted.values

  if (
    !isTRUE(fit$converged) ||
      fit$rank < n_parameters ||
      anyNA(coefficients) ||
      any(!is.finite(coefficients)) ||
      any(!is.finite(fitted_mean) | fitted_mean <= 0)
  ) {
    stop("The capture-model M-step failed to converge.", call. = FALSE)
  }

  structure(
    list(
      coefficients = coefficients,
      fitted_mean = fitted_mean,
      objective = sum(cell_counts * log(fitted_mean) - fitted_mean),
      iterations = fit$iter,
      converged = TRUE
    ),
    class = "crc_capture_mstep"
  )

}