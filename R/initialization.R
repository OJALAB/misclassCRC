#' Initialize Latent-State Weights and Capture-Cell Counts
#'
#' @import data.table
#' @importFrom stats rgamma
#' 
#' @noRd
initialize_crc <- function(
  model_matrices,
  misclass,
  latent_classes,
  control
) {

  init_alpha <- if (is.null(control$init_alpha)) 20 else control$init_alpha

  states <- model_matrices$states
  n_states <- nrow(states$data)

  if (latent_classes == 1L) {
    latent_weight <- rep(1, n_states)
  } else {
    misclass_names <- if (is.null(misclass)) character(0L) else misclass$variables
    groups <- copy(states$data[, misclass_names, with = FALSE])
    groups[, c(".observation_id", ".state_row") :=
      list(states$observation_id, seq_len(n_states))]
    group_names <- c(".observation_id", misclass_names)

    allocations <- groups[, {
      draw <- rgamma(.N, shape = init_alpha)
      list(.state_row = get(".state_row"),
           .latent_weight = draw / sum(draw))
    }, by = group_names]
    setorderv(allocations, ".state_row")
    latent_weight <- allocations[[".latent_weight"]]
  }

  weights <- data.table(
    .observation_id = states$observation_id,
    .raw = states$misclass_probability * latent_weight
  )
  normalizers <- weights[, sum(get(".raw")), by = ".observation_id"][["V1"]]
  if (any(!is.finite(normalizers) | normalizers <= 0)) {
    stop("Some observations have no admissible initial states.", call. = FALSE)
  }
  weights[, ".state_weight" :=
    get(".raw") / sum(get(".raw")), by = ".observation_id"]

  n_cells <- nrow(model_matrices$capture$cells)

  cell_counts <- numeric(n_cells)

  aggregation_data <- data.table(
    .cell = model_matrices$capture_index,
    .weight = weights[[".state_weight"]]
  )

  aggregated <- aggregation_data[,
    list(count = sum(get(".weight"))),
    by = ".cell"
  ]

  cell_counts[aggregated[[".cell"]]] <- aggregated[["count"]]

  cell_counts[model_matrices$capture$unobserved] <- 1

  structure(
    list(state_weights = weights[[".state_weight"]],
         cell_counts = cell_counts, init_alpha = init_alpha),
    class = "crc_initialization"
  )

}
