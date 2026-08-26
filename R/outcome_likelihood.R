#' Subtract Probabilities on the Log Scale
#' 
#' @noRd
log_sub_crc <- function(log_x, log_y) {

  difference <- log_y - log_x

  both_zero <- is.infinite(log_x) &
    log_x < 0 &
    is.infinite(log_y) &
    log_y < 0

  difference[both_zero] <- -Inf

  log_x + log(-expm1(difference))

}

#' Calculate Outcome Log-Likelihood Contributions
#' 
#' @importFrom stats dnbinom pnbinom dpois ppois
#' 
#' @noRd
outcome_loglik_crc <- function(
  outcome,
  mu,
  outcome_dist,
  sigma = NULL
) {

  n <- nrow(outcome)

  if (
    length(mu) != n ||
      anyNA(mu) ||
      any(!is.finite(mu)) ||
      any(mu <= 0)
  ) {
    stop(
     "`mu` must contain one finite positive value per outcome state.",
    call. = FALSE
    )
  }

  if (outcome_dist == "ztnegbin") {
    if (
      is.null(sigma) ||
        !is.numeric(sigma) ||
        length(sigma) != 1L ||
        anyNA(sigma) ||
        any(!is.finite(sigma)) ||
        any(sigma <= 0)
    ) {
      stop("`sigma` must be a finite positive number.", call. = FALSE)
    }
    
    size <- 1 / sigma
    log_mass <- function(q) dnbinom(q, mu = mu, size = size, log = TRUE)
    log_survival <- function(q) {
      pnbinom(q, mu = mu, size = size, lower.tail = FALSE, log.p = TRUE)
    }
  } else if (outcome_dist == "ztpois") {
    log_mass <- function(q) dpois(q, lambda = mu, log = TRUE)
    log_survival <- function(q) {
      ppois(q, lambda = mu, lower.tail = FALSE, log.p = TRUE)
    }
  }

  log_normalizer <- log_survival(0)
  loglik <- numeric(n)

  exact <- outcome$status == "exact"
  interval <- outcome$status == "interval_censored"
  right <- outcome$status == "right_censored"

  loglik[exact] <- log_mass(outcome$lower)[exact] - log_normalizer[exact]

  loglik[interval] <- log_sub_crc(
    log_survival(outcome$lower - 1)[interval],
    log_survival(outcome$upper)[interval]
  ) - log_normalizer[interval]

  loglik[right] <- log_survival(outcome$lower - 1)[right] -
    log_normalizer[right]

  loglik

}