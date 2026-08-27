#' Perform the E-Step
#' 
#' @import data.table
#' 
#' @noRd
expectation_crc <- function(
  model_matrices,
  capture_mean,
  outcome_loglik = NULL
) {

  states <- model_matrices$states
  n_states <- nrow(states$data)
  n_cells <- nrow(model_matrices$capture$cells)

  if (length(capture_mean) != n_cells) {
    stop(
      "Capture-model fitted means are not aligned with the capture cells.",
      call. = FALSE
    )
  }

  if (
    anyNA(capture_mean) ||
      any(!is.finite(capture_mean)) ||
      any(capture_mean <= 0)
  ) {
    stop(
      "The capture model produced non-finite or non-positive fitted means.",
      call. = FALSE
    )
  }

  if (is.null(outcome_loglik)) {
    outcome_loglik <- numeric(n_states)
  }

  if (
    length(outcome_loglik) != n_states ||
      anyNA(outcome_loglik) ||
      any(outcome_loglik == Inf)
  ) {
    stop(
      "The outcome model produced invalid log-likelihood contributions.",
      call. = FALSE
    )
  }

  log_misclass <- rep(-Inf, n_states)
  positive <- states$misclass_probability > 0
  log_misclass[positive] <- log(states$misclass_probability[positive])

  state_weights <- data.table(
    .observation_id = states$observation_id,
    .log_weight = log(capture_mean[model_matrices$capture_index]) +
      log_misclass +
      outcome_loglik
  )

  admissible <- state_weights[,
    any(is.finite(get(".log_weight"))),
    by = ".observation_id"
  ][["V1"]]
  if (!all(admissible)) {
    stop("Some observations have no admissible states.", call. = FALSE)
  }

  state_weights[,
    ".state_weight" := {
      z <- get(".log_weight")
      z_max <- max(z)
      scaled <- exp(z - z_max)
      scaled / sum(scaled)
    },
    by = ".observation_id"
  ]

  cell_counts <- numeric(n_cells)
  aggregated <- data.table(
    .cell = model_matrices$capture_index,
    .weight = state_weights[[".state_weight"]]
  )[, list(.count = sum(get(".weight"))), by = ".cell"]

  cell_counts[aggregated[[".cell"]]] <- aggregated[[".count"]]
  unobserved <- model_matrices$capture$unobserved
  cell_counts[unobserved] <- capture_mean[unobserved]

  structure(
    list(
      state_weights = state_weights[[".state_weight"]],
      cell_counts = cell_counts
    ),
    class = "crc_estep"
  )

}