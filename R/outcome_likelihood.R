#' Subtract Probabilities on the Log Scale
#' 
#' @noRd
log_sub_crc <- function(
  log_x,
  log_y
) {

  difference <- log_y - log_x

  both_zero <- is.infinite(log_x) &
    log_x < 0 &
    is.infinite(log_y) &
    log_y < 0

  difference[both_zero] <- -Inf

  log_x + log(-expm1(difference))

}

#' Calculate log(1 - exp(x))
#'
#' @noRd
log1mexp_crc <- function(
  x
) {
  out <- numeric(length(x))
  left <- x < -log(2)
  out[left] <- log1p(-exp(x[left]))
  out[!left] <- log(-expm1(x[!left]))
  out
}

#' Check Whether to Use the Poisson Limit
#' 
#' @noRd
nb_use_poisson_crc <- function(
  mu,
  sigma
) {
  sigma * mu <= sqrt(.Machine$double.eps)
}

#' Calculate NB Probabilities on the Log Scale
#'
#' @noRd
nb_log_mass_crc <- function(
  q,
  mu,
  sigma
) {
  q <- rep_len(q, length(mu))
  poisson <- nb_use_poisson_crc(mu, sigma)
  out <- numeric(length(mu))

  out[poisson] <- dpois(q[poisson], mu[poisson], log = TRUE)
  out[!poisson] <- dnbinom(
    q[!poisson],
    mu = mu[!poisson],
    size = 1 / sigma,
    log = TRUE
  )
  out
}

#' Calculate NB Upper-Tail Probabilities on the Log Scale
#'
#' @noRd
nb_log_survival_crc <- function(
  q,
  mu,
  sigma
) {
  q <- rep_len(q, length(mu))
  out <- rep(NA_real_, length(mu))
  valid <- !is.na(q)
  poisson <- valid & nb_use_poisson_crc(mu, sigma)

  out[poisson] <- ppois(
    q[poisson],
    mu[poisson],
    lower.tail = FALSE,
    log.p = TRUE
  )

  zero <- valid & !poisson & q == 0
  log_p0 <- -log1p(sigma * mu[zero]) / sigma
  out[zero] <- log1mexp_crc(log_p0)

  ordinary <- valid & !poisson & !zero
  out[ordinary] <- pnbinom(
    q[ordinary],
    mu = mu[ordinary],
    size = 1 / sigma,
    lower.tail = FALSE,
    log.p = TRUE
  )
  out
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
    
    log_mass <- function(q) nb_log_mass_crc(q, mu, sigma)
    log_survival <- function(q) nb_log_survival_crc(q, mu, sigma)
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