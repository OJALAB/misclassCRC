#' Build the Capture-Model Matrix
#' 
#' @import data.table
#' @importFrom stats model.matrix setNames
#'
#' @noRd
build_capture_matrix <- function(
  data,
  captures,
  capture_model,
  misclass,
  latent_classes
) {
  capture_names <- captures$names
  categorical_names <- setdiff(
    capture_model$variables,
    capture_names
  )

  get_levels <- function(variable) {
    if (
      !is.null(misclass) &&
        variable %in% misclass$variables
    ) {
      return(misclass$specifications[[variable]]$levels)
    }

    x <- data[[variable]]

    if (is.factor(x)) levels(x) else unique(x)
  }

  grid_values <- setNames(
    rep(list(0:1), length(capture_names)),
    capture_names
  )

  for (variable in categorical_names) {
    levels <- get_levels(variable)
    grid_values[[variable]] <- factor(levels, levels = levels)
  }

  grid_values[[".latent"]] <- factor(seq_len(latent_classes))

  cells <- do.call(
    CJ,
    c(grid_values, list(sorted = FALSE, unique = TRUE))
  )

  matrix <- model.matrix(
    capture_model$terms,
    data = cells
  )

  list(
    cells = cells,
    matrix = matrix,
    unobserved = rowSums(
      as.matrix(cells[, capture_names, with = FALSE])
    ) == 0L
  )

}

#' Build Initial Observation States
#' 
#' @import data.table
#' 
#' @noRd
build_observation_states <- function(
  data,
  captures,
  capture_model,
  outcome,
  outcome_model,
  misclass,
  latent_classes
) {

  outcome_names <- if (is.null(outcome_model)) {
    character(0L)
  } else {
    outcome_model$variables
  }

  misclass_names <- if (is.null(misclass)) {
    character(0L)
  } else {
    misclass$variables
  }

  variable_names <- unique(c(
    captures$names,
    capture_model$variables,
    outcome_names
  ))

  observation_id <- rep(
    seq_len(nrow(data)),
    each = latent_classes
  )

  state_data <- as.data.table(data)[
    observation_id,
    variable_names,
    with = FALSE
  ]

  state_data[,
    (".latent") := factor(
      rep(seq_len(latent_classes), times = nrow(data)),
      levels = seq_len(latent_classes)
    )
  ]

  observed_misclass <- if (length(misclass_names) == 0L) {
    NULL
  } else {
    as.data.table(data)[
      observation_id,
      misclass_names,
      with = FALSE
    ]
  }

  structure(
    list(
      data = state_data,
      observation_id = observation_id,
      observed_misclass = observed_misclass,
      outcome = if (is.null(outcome)) NULL else outcome[observation_id]
    ),
    class = "crc_observation_states"
  )

}

#' Expand States for a Misclassified Variable
#' 
#' @import data.table
#' 
#' @noRd
expand_misclass_states <- function(
  states,
  variable,
  misclass
) {

  specification <- misclass$specifications[[variable]]
  true_levels <- specification$levels
  n_levels <- length(true_levels)

  row_index <- rep(
    seq_len(nrow(states$data)),
    each = n_levels
  )

  candidate <- rep(
    true_levels,
    times = nrow(states$data)
  )

  observed <- as.character(
    states$observed_misclass[[variable]][row_index]
  )

  probability <- specification$matrix[cbind(
    match(candidate, rownames(specification$matrix)),
    match(observed, colnames(specification$matrix))
  )]

  known_truth <- specification$is_true[
    states$observation_id[row_index]
  ]

  probability[known_truth] <- as.numeric(
    candidate[known_truth] == observed[known_truth]
  )

  previous_probability <- states$misclass_probability

  if (is.null(previous_probability)) {
    previous_probability <- rep(1, nrow(states$data))
  }

  state_data <- states$data[row_index]
  state_data[,
    (variable) := factor(
      candidate,
      levels = true_levels
    )
  ]

  structure(
    list(
      data = state_data,
      observation_id = states$observation_id[row_index],
      observed_misclass = states$observed_misclass[row_index],
      outcome = if (is.null(states$outcome)) {
        NULL
      } else {
        states$outcome[row_index]
      },
      misclass_probability = previous_probability[row_index] * probability
    ),
    class = "crc_observation_states"
  )

}

#' Expand All Misclassified Variables
#'
#' @noRd
expand_all_misclass_states <- function(
  states,
  misclass
) {

  if (is.null(misclass)) {
    states$misclass_probability <- rep(
      1,
      nrow(states$data)
    )

    return(states)
  }

  for (variable in misclass$variables) {
    states <- expand_misclass_states(
      states = states,
      variable = variable,
      misclass = misclass
    )
  }

  states
  
}

#' Build the Outcome-Model Matrix
#'
#' @importFrom stats model.matrix
#'
#' @noRd
build_outcome_matrix <- function(
  states,
  outcome_model
) {

  if (is.null(outcome_model)) {
    return(NULL)
  }

  matrix <- model.matrix(
    outcome_model$terms,
    data = states$data
  )

  n_states <- nrow(states$data)

  aligned <- (nrow(matrix) == n_states &&
    length(states$observation_id) == n_states &&
    !is.null(states$outcome) &&
    nrow(states$outcome) == n_states &&
    length(states$misclass_probability) == n_states)

  if (!aligned) {
    stop(
      "Outcome-model components are not aligned with observation states.",
      call. = FALSE
    )
  }

  matrix

}

#' Build Capture-Cell Indices
#' 
#' @noRd
build_capture_index <- function(
  states,
  capture,
  captures,
  capture_model
) {

  key_names <- unique(c(
    captures$names,
    capture_model$variables,
    ".latent"
  ))

  lookup <- copy(
    capture$cells[, key_names, with = FALSE]
  )
  lookup[
    ,
    ("capture_index") := seq_len(nrow(lookup))
  ]

  index <- lookup[
    states$data,
    on = key_names
  ][["capture_index"]]

  if (anyNA(index)) {
    stop(
      "Some observation states do not match the capture-cell grid.",
      call. = FALSE
    )
  }

  index

}

#' Build Model Matrices
#' 
#' @noRd
build_model_matrices <- function(
  data,
  captures,
  capture_model,
  outcome,
  outcome_model,
  misclass,
  latent_classes
) {

  capture <- build_capture_matrix(
    data,
    captures,
    capture_model,
    misclass,
    latent_classes
  )

  states <- build_observation_states(
    data,
    captures,
    capture_model,
    outcome,
    outcome_model,
    misclass,
    latent_classes
  )
  states <- expand_all_misclass_states(states, misclass)

  structure(
    list(
      capture = capture,
      states = states,
      outcome_matrix = build_outcome_matrix(states, outcome_model),
      capture_index = build_capture_index(
        states,
        capture,
        captures,
        capture_model
      )
    ),
    class = "crc_model_matrices"
  )

}