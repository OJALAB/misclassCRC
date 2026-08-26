library(misclassCRC)

outcome <- data.frame(
  lower = c(2, 3, 4, NA),
  upper = c(2, 5, Inf, NA),
  status = c("exact", "interval_censored", "right_censored", "missing")
)
mu <- c(1.5, 3, 6, 2)   
pois_normalizer <- ppois(0, mu, lower.tail = FALSE, log.p = TRUE)
expected_pois <- c(
  dpois(2, mu[1L], log = TRUE) - pois_normalizer[1L],
  log(ppois(5, mu[2L]) - ppois(2, mu[2L])) - pois_normalizer[2L],
  ppois(3, mu[3L], lower.tail = FALSE, log.p = TRUE) - pois_normalizer[3L],
  0
)
result_pois <- misclassCRC:::outcome_loglik_crc(outcome, mu, "ztpois")
expect_equal(result_pois, expected_pois)

sigma <- 0.4
size <- 1 / sigma
nb_normalizer <- pnbinom(0, mu = mu, size = size, lower.tail = FALSE,
  log.p = TRUE)
expected_nb <- c(
  dnbinom(2, mu = mu[1L], size = size, log = TRUE) - nb_normalizer[1L],
  log(pnbinom(5, mu = mu[2L], size = size) -
    pnbinom(2, mu = mu[2L], size = size)) - nb_normalizer[2L],
  pnbinom(3, mu = mu[3L], size = size, lower.tail = FALSE,
    log.p = TRUE) - nb_normalizer[3L],
  0
)
result_nb <- misclassCRC:::outcome_loglik_crc(outcome, mu, "ztnegbin", sigma)
expect_equal(result_nb, expected_nb)

far_exact <- data.frame(lower = 500, upper = 500, status = "exact")
far_interval <- data.frame(lower = 500, upper = 501,
  status = "interval_censored")
expect_true(is.finite(misclassCRC:::outcome_loglik_crc(far_exact, 10, "ztpois")))
expect_true(is.finite(misclassCRC:::outcome_loglik_crc(far_interval, 10, "ztpois")))
right_from_one <- data.frame(lower = 1, upper = Inf,
  status = "right_censored")
expect_equal(misclassCRC:::outcome_loglik_crc(right_from_one, 2, "ztpois"), 0)

expect_error(misclassCRC:::outcome_loglik_crc(outcome, mu[-1L], "ztpois"),
  pattern = "one finite positive value")
expect_error(misclassCRC:::outcome_loglik_crc(outcome, c(0, mu[-1L]), "ztpois"),
  pattern = "one finite positive value")
expect_error(misclassCRC:::outcome_loglik_crc(outcome, mu, "ztnegbin"),
  pattern = "finite positive number")
expect_error(misclassCRC:::outcome_loglik_crc(outcome, mu, "ztnegbin", c(1, 2)),
  pattern = "finite positive number")
