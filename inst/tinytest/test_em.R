library(misclassCRC)

data <- data.frame(
  source_1 = rep(c(1L, 0L, 0L, 1L, 1L, 0L), 4L),
  source_2 = rep(c(0L, 1L, 0L, 1L, 0L, 1L), 4L),
  source_3 = rep(c(0L, 0L, 1L, 0L, 1L, 1L), 4L),
  group = factor(rep(c("a", "b"), 12L)),
  outcome = rep(c(1, 2, 3, 10, 2, 14), 4L)
)

prepare_case <- function(outcome_dist, with_outcome = TRUE, control = NULL) {
  misclassCRC:::prepare_crc(
    data, ~ source_1 + source_2 + source_3,
    ~ source_1 + source_2 + source_3 + group,
    outcome = if (with_outcome) "outcome" else NULL,
    outcome_formula = if (with_outcome) ~ group else NULL,
    outcome_dist = outcome_dist, misclass = NULL, latent_classes = 1L,
    control = control, verbose = FALSE
  )
}

preparation <- prepare_case("ztnegbin")
first <- misclassCRC:::em_iteration_crc(
  preparation$model_matrices, preparation$initialization$cell_counts,
  preparation$initialization$state_weights, preparation$outcome_dist,
  preparation$control, NULL, NULL
)
expect_inherits(first, "crc_em_iteration")
expect_inherits(first$capture_fit, "crc_capture_mstep")
expect_inherits(first$outcome_fit, "crc_outcome_mstep")
expect_true(is.finite(first$outcome_fit$sigma) && first$outcome_fit$sigma > 0)
expect_equal(as.numeric(tapply(first$state_weights,
  preparation$model_matrices$states$observation_id, sum)), rep(1, nrow(data)))

outcome_start <- c(first$outcome_fit$coefficients,
  log_sigma = log(first$outcome_fit$sigma))
second <- misclassCRC:::em_iteration_crc(
  preparation$model_matrices, first$cell_counts, first$state_weights,
  preparation$outcome_dist, preparation$control,
  first$capture_fit$coefficients, outcome_start
)
expect_true(all(is.finite(second$capture_fit$fitted_mean)))
expect_true(all(is.finite(second$outcome_fit$fitted_mean)))
expect_false(second$outcome_fit$retried)

fit <- misclassCRC:::perform_em_crc(
  preparation$model_matrices, preparation$initialization,
  preparation$outcome_dist, preparation$control, FALSE
)
expect_inherits(fit, "crc_em_fit")
expect_true(fit$converged)
expect_true(fit$iterations > 1L)
expect_true(fit$convergence$maximum < preparation$control$em$tolerance)
expect_equal(names(fit$convergence$changes),
  c("state_weights", "cell_counts"))
expect_true(all(is.finite(fit$convergence$changes)))

poisson <- prepare_case("ztpois")
poisson_fit <- misclassCRC:::perform_em_crc(
  poisson$model_matrices, poisson$initialization,
  poisson$outcome_dist, poisson$control, FALSE
)
expect_true(poisson_fit$converged)
expect_null(poisson_fit$outcome_fit$sigma)
expect_equal(names(poisson_fit$convergence$changes),
  c("state_weights", "cell_counts"))
expect_true(all(is.finite(poisson_fit$convergence$changes)))

without_outcome <- prepare_case("ztpois", with_outcome = FALSE)
without_outcome_fit <- misclassCRC:::perform_em_crc(
  without_outcome$model_matrices, without_outcome$initialization,
  without_outcome$outcome_dist, without_outcome$control, FALSE
)
expect_true(without_outcome_fit$converged)
expect_null(without_outcome_fit$outcome_fit)
expect_equal(names(without_outcome_fit$convergence$changes),
  c("state_weights", "cell_counts"))
expect_true(all(is.finite(without_outcome_fit$convergence$changes)))

limited <- prepare_case("ztnegbin", control = list(em = list(max_iter = 1L)))
limited_fit <- misclassCRC:::perform_em_crc(
  limited$model_matrices, limited$initialization,
  limited$outcome_dist, limited$control, FALSE
)
expect_false(limited_fit$converged)
expect_equal(limited_fit$iterations, 1L)
expect_equal(limited_fit$convergence$maximum, Inf)
