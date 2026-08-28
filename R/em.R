#' Perform One EM Iteration
#' 
#' @noRd
em_iteration_crc <- function(
  model_matrices,
  cell_counts,
  state_weights,
  outcome_dist,
  control,
  capture_start,
  outcome_start
) {

  capture_fit <- maximize_capture_crc(
    model_matrices = model_matrices,
    cell_counts = cell_counts,
    start = capture_start,
    glm_control = list(
      maxit = control$capture$max_iter,
      epsilon = control$capture$tolerance
    )
  )

  outcome_fit <- NULL
  outcome_loglik <- NULL

  if (!is.null(model_matrices$outcome_matrix)) {
    outcome_fit <- maximize_outcome_crc(
      model_matrices = model_matrices,
      state_weights = state_weights,
      outcome_dist = outcome_dist,
      start = outcome_start,
      optim_control = list(
        maxit = control$outcome$max_iter,
        reltol = control$outcome$relative_tolerance
      )
    )

    outcome_loglik <- outcome_loglik_crc(
      outcome = model_matrices$states$outcome,
      mu = outcome_fit$fitted_mean,
      outcome_dist = outcome_dist,
      sigma = outcome_fit$sigma
    )
  }

  expectation <- expectation_crc(
    model_matrices = model_matrices,
    capture_mean = capture_fit$fitted_mean,
    outcome_loglik = outcome_loglik
  )

  structure(
    list(
      capture_fit = capture_fit,
      outcome_fit = outcome_fit,
      state_weights = expectation$state_weights,
      cell_counts = expectation$cell_counts
    ),
    class = "crc_em_iteration"
  )

}

#' Perform the EM Algorithm
#' 
#' @noRd
perform_em_crc <- function(
  model_matrices,
  initialization,
  outcome_dist,
  control,
  verbose
) {

  state_weights <- initialization$state_weights
  cell_counts <- initialization$cell_counts
  previous <- NULL
  converged <- FALSE

  for (iter in seq_len(control$em$max_iter)) {
    capture_start <- if (is.null(previous)) {
      NULL
    } else {
      previous$capture_fit$coefficients
    }

    outcome_start <- NULL
    if (!is.null(previous$outcome_fit)) {
      outcome_start <- previous$outcome_fit$coefficients
      if (outcome_dist == "ztnegbin") {
        outcome_start <- c(
          outcome_start,
          log_sigma = log(previous$outcome_fit$sigma)
        )
      }
    }

    current <- em_iteration_crc(
      model_matrices = model_matrices,
      cell_counts = cell_counts,
      state_weights = state_weights,
      outcome_dist = outcome_dist,
      control = control,
      capture_start = capture_start,
      outcome_start = outcome_start
    )

    changes <- c(
      state_weights = Inf,
      capture_mean = Inf,
      outcome_mean = if (is.null(current$outcome_fit)) NA_real_ else Inf,
      sigma = if (is.null(current$outcome_fit$sigma)) NA_real_ else Inf
    )

    if (!is.null(previous)) {
      changes["state_weights"] <- max(abs(
        current$state_weights - previous$state_weights
      ))
      changes["capture_mean"] <- max(abs(
        log(current$capture_fit$fitted_mean) -
          log(previous$capture_fit$fitted_mean)
      ))
      if (!is.null(current$outcome_fit)) {
        changes["outcome_mean"] <- max(abs(
          log(current$outcome_fit$fitted_mean) -
            log(previous$outcome_fit$fitted_mean)
        ))
      }
      if (!is.null(current$outcome_fit$sigma)) {
        changes["sigma"] <- abs(
          log(current$outcome_fit$sigma) -
            log(previous$outcome_fit$sigma)
        )
      }
    }

    maximum <- max(changes, na.rm = TRUE)
    converged <- is.finite(maximum) &&
      maximum < control$em$tolerance

    previous <- current
    state_weights <- current$state_weights
    cell_counts <- current$cell_counts

    if (converged) {
      break
    }
  }

  current$iterations <- iter
  current$converged <- converged
  current$convergence <- list(
    maximum = maximum,
    changes = changes
  )
  class(current) <- "crc_em_fit"
  current

}