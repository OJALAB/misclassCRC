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

#' Maximize the Outcome-Model Objective
#'
#' @importFrom stats glm.fit optim poisson
#'
#' @noRd
maximize_outcome_crc <- function(
  model_matrices,
  state_weights,
  outcome_dist,
  start = NULL,
  optim_control = list()
) {
  matrix <- model_matrices$outcome_matrix
  outcome <- model_matrices$states$outcome
  observation_id <- model_matrices$states$observation_id

  if (
    is.null(matrix) ||
      is.null(outcome) ||
      nrow(outcome) != nrow(matrix) ||
      length(observation_id) != nrow(matrix) ||
      nrow(matrix) != length(state_weights) ||
      anyNA(state_weights) ||
      any(!is.finite(state_weights)) ||
      any(state_weights < 0)
  ) {
    stop("Outcome-model inputs are invalid.", call. = FALSE)
  }

  n_parameters <- ncol(matrix)
  active <- state_weights > 0
  informative <- outcome$status != "missing" & active

  if (!any(informative)) {
    stop("No informative outcomes have positive state weight.", call. = FALSE)
  }
  if (qr(matrix[informative, , drop = FALSE])$rank < n_parameters) {
    stop("The outcome-model matrix is rank deficient.", call. = FALSE)
  }

  exact <- outcome$status == "exact" & active
  n_exact <- length(unique(observation_id[exact]))
  beta_start <- numeric(n_parameters)
  names(beta_start) <- colnames(matrix)
  initial_fit <- NULL

  if (n_exact >= n_parameters) {
    initial_fit <- tryCatch(
      suppressWarnings(glm.fit(
        x = matrix[exact, , drop = FALSE],
        y = outcome$lower[exact],
        weights = state_weights[exact],
        family = poisson("log")
      )),
      error = function(error) NULL
    )
  }

  valid_initial_fit <- !is.null(initial_fit) &&
    isTRUE(initial_fit$converged) &&
    initial_fit$rank == n_parameters &&
    all(is.finite(initial_fit$coefficients))

  if (valid_initial_fit) {
    beta_start <- initial_fit$coefficients
  } else {
    w <- state_weights[informative]
    initial_mean <- sum(w * outcome$lower[informative]) / sum(w)
    beta_start[1L] <- log(initial_mean)
  }

  default_start <- beta_start
  if (outcome_dist == "ztnegbin") {
    sigma_start <- 0.25
    if (n_exact >= 2L) {
      w <- state_weights[exact]
      mean_y <- sum(w * outcome$lower[exact]) / sum(w)
      variance_y <- sum(w * (outcome$lower[exact] - mean_y)^2) / sum(w)
      sigma_start <- max((variance_y - mean_y) / mean_y^2, 1e-4)
    }
    default_start <- c(default_start, log_sigma = log(sigma_start))
  }
  if (is.null(start)) {
    start <- default_start
  }
  if (
    length(start) != length(default_start) ||
      any(!is.finite(start))
  ) {
    stop("Outcome-model starting values are invalid.", call. = FALSE)
  }

  expected_names <- names(default_start)

  if (is.null(names(start))) {
    names(start) <- expected_names
  } else {
    if (
      anyDuplicated(names(start)) || !setequal(names(start), expected_names)
    ) {
      stop("Outcome-model starting-value names are invalid.", call. = FALSE)
    }
    start <- start[expected_names]
  }

  penalty <- 1e100

  objective <- function(parameters) {
    beta <- parameters[seq_len(n_parameters)]
    mu <- exp(drop(matrix %*% beta))
    sigma <- if (outcome_dist == "ztnegbin") {
      exp(parameters[n_parameters + 1L])
    } else {
      NULL
    }

    loglik <- tryCatch(
      outcome_loglik_crc(outcome, mu, outcome_dist, sigma),
      error = function(error) NULL
    )
    if (is.null(loglik) || any(!is.finite(loglik[active]))) {
      return(penalty)
    }

    value <- -sum(state_weights[active] * loglik[active])
    if (is.finite(value)) value else penalty
  }

  run_optim <- function(candidate) {
    initial <- objective(candidate)
    if (initial >= penalty) {
      return(NULL)
    }

    result <- tryCatch(
      optim(candidate, objective, method = "BFGS", control = optim_control),
      error = function(error) NULL
    )
    tolerance <- 1e-8 * (1 + abs(initial))
    if (
      is.null(result) ||
        result$convergence != 0L ||
        !is.finite(result$value) ||
        result$value >= penalty ||
        result$value > initial + tolerance
    ) {
      return(NULL)
    }
    result
  }

  fit <- run_optim(start)
  retried <- FALSE

  if (
    is.null(fit) &&
      !isTRUE(all.equal(unname(start), unname(default_start)))
  ) {
    fit <- run_optim(default_start)
    retried <- TRUE
  }

  if (is.null(fit)) {
    stop("The outcome-model M-step failed to converge.", call. = FALSE)
  }

  beta <- fit$par[seq_len(n_parameters)]
  fitted_mean <- exp(drop(matrix %*% beta))
  sigma <- if (outcome_dist == "ztnegbin") {
    unname(exp(fit$par[n_parameters + 1L]))
  } else {
    NULL
  }

  if (
    any(!is.finite(beta)) ||
      any(!is.finite(fitted_mean) | fitted_mean <= 0) ||
      (!is.null(sigma) && (!is.finite(sigma) || sigma <= 0))
  ) {
    stop("The outcome-model M-step produced invalid estimates.", call. = FALSE)
  }

  structure(
    list(
      coefficients = beta,
      fitted_mean = fitted_mean,
      sigma = if (outcome_dist == "ztnegbin") {
        sigma
      } else {
        NULL
      },
      objective = -fit$value,
      evaluations = fit$counts,
      retried = retried,
      converged = TRUE
    ),
    class = "crc_outcome_mstep"
  )

}
