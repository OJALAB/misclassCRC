#' @noRd
validate_crc_input <- function(
  data,
  captures,
  capture_formula,
  outcome,
  outcome_formula,
  outcome_dist,
  misclass,
  latent_classes,
  control
) {

  if (
    !is.data.frame(data)
  ) {
    stop("`data` must be a `data.frame` or `data.table`", call. = FALSE)
  }

  if (
    !inherits(captures, "formula") ||
      length(captures) != 2L
  ) {
    stop("`captures` must be a one-sided formula.", call. = FALSE)
  }

  if (
    !inherits(capture_formula, "formula") ||
      length(capture_formula) != 2L
  ) {
    stop("`capture_formula` must be a one-sided formula.", call. = FALSE)
  }

  if (is.null(outcome) && !is.null(outcome_formula)) {
    stop(
      "`outcome_formula` must be `NULL` when `outcome` is `NULL`.",
      call. = FALSE
    )
  }

  if (
    !is.null(outcome) &&
      (is.null(outcome_formula) ||
        !inherits(outcome_formula, "formula") ||
        length(outcome_formula) != 2L)
  ) {
    stop(
      "`outcome_formula` must be supplied as a one-sided formula.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(latent_classes) ||
      length(latent_classes) != 1L ||
      is.na(latent_classes) ||
      !is.finite(latent_classes) ||
      latent_classes < 1 ||
      latent_classes != floor(latent_classes)
  ) {
    stop(
      "`latent_classes` must be a positive integer.",
      call. = FALSE
    )
  }

  invisible(TRUE)

}

#' @noRd
validate_categorical_variables <- function(
  variables,
  data,
  context
) {

  variables_with_missing <- variables[
    vapply(
      variables,
      function(variable) anyNA(data[[variable]]),
      logical(1L)
    )
  ]

  if (length(variables_with_missing) > 0L) {
    stop(
      "The following variables in `",
      context,
      "` contain missing values: ",
      paste0("`", variables_with_missing, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  non_categorical_variables <- variables[
    !vapply(
      variables,
      function(variable) {
        is.factor(data[[variable]]) ||
          is.character(data[[variable]])
      },
      logical(1L)
    )
  ]

  if (length(non_categorical_variables) > 0L) {
    stop(
      "The following variables in `",
      context,
      "` must be factor or character variables: ",
      paste0("`", non_categorical_variables, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
  
}