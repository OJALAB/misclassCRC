library(misclassCRC)

data <- data.frame(
  source_1 = c(1L, 0L, 0L, 1L),
  source_2 = c(0L, 1L, 1L, 0L),
  source_3 = c(0L, 0L, 1L, 1L)
)
fit <- crc_fit(data, ~ source_1 + source_2 + source_3,
  ~ source_1 + source_2 + source_3)
matrices <- fit$model_matrices
counts <- fit$initialization$cell_counts
result <- misclassCRC:::maximize_capture_crc(matrices, counts)

expect_inherits(result, "crc_capture_mstep")
expect_true(result$converged)
expect_equal(length(result$coefficients), ncol(matrices$capture$matrix))
expect_equal(length(result$fitted_mean), length(counts))
expect_true(all(is.finite(result$coefficients)))
expect_true(all(is.finite(result$fitted_mean) & result$fitted_mean > 0))
expect_equal(result$objective,
  sum(counts * log(result$fitted_mean) - result$fitted_mean))

warm <- misclassCRC:::maximize_capture_crc(
  matrices, counts, start = result$coefficients)
expect_equal(warm$fitted_mean, result$fitted_mean)
expect_error(misclassCRC:::maximize_capture_crc(matrices, counts[-1L]),
  pattern = "cell counts")
expect_error(misclassCRC:::maximize_capture_crc(matrices, replace(counts, 1L, -1)),
  pattern = "cell counts")
expect_error(misclassCRC:::maximize_capture_crc(matrices, counts, start = 0),
  pattern = "starting values")

rank_deficient <- matrices
rank_deficient$capture$matrix <- cbind(matrices$capture$matrix,
  duplicate = matrices$capture$matrix[, 1L])
expect_error(misclassCRC:::maximize_capture_crc(rank_deficient, counts),
  pattern = "failed to converge")

outcome_data <- data.frame(
  source_1 = rep(c(1L, 0L, 0L, 1L, 1L, 0L), 4L),
  source_2 = rep(c(0L, 1L, 0L, 1L, 0L, 1L), 4L),
  source_3 = rep(c(0L, 0L, 1L, 0L, 1L, 1L), 4L),
  group = factor(rep(c("a", "b"), 12L)),
  outcome = rep(c(1, 2, 3, 10, 2, 14), 4L)
)
outcome_fit <- crc_fit(
  outcome_data, ~ source_1 + source_2 + source_3,
  ~ source_1 + source_2 + source_3 + group,
  outcome = "outcome", outcome_formula = ~ group
)
outcome_matrices <- outcome_fit$model_matrices
state_weights <- outcome_fit$initialization$state_weights

poisson_result <- misclassCRC:::maximize_outcome_crc(
  outcome_matrices, state_weights, "ztpois")
negbin_result <- misclassCRC:::maximize_outcome_crc(
  outcome_matrices, state_weights, "ztnegbin")
expect_inherits(poisson_result, "crc_outcome_mstep")
expect_true(poisson_result$converged)
expect_null(poisson_result$sigma)
expect_equal(names(poisson_result$coefficients),
  colnames(outcome_matrices$outcome_matrix))
expect_true(all(is.finite(poisson_result$fitted_mean) &
  poisson_result$fitted_mean > 0))
expect_true(negbin_result$converged)
expect_true(is.finite(negbin_result$sigma) && negbin_result$sigma > 0)
expect_null(names(negbin_result$sigma))
expect_true(is.finite(negbin_result$objective))

warm_start <- c(log_sigma = log(negbin_result$sigma),
  rev(negbin_result$coefficients))
warm_result <- misclassCRC:::maximize_outcome_crc(
  outcome_matrices, state_weights, "ztnegbin", start = warm_start)
expect_false(warm_result$retried)
expect_equal(warm_result$objective, negbin_result$objective, tolerance = 1e-6)

bad_start <- setNames(rep(1000, 3L),
  c(colnames(outcome_matrices$outcome_matrix), "log_sigma"))
retry_result <- misclassCRC:::maximize_outcome_crc(
  outcome_matrices, state_weights, "ztnegbin", start = bad_start)
expect_true(retry_result$retried)
expect_error(misclassCRC:::maximize_outcome_crc(
  outcome_matrices, state_weights[-1L], "ztpois"), pattern = "inputs")

rank_data <- outcome_data
rank_data$outcome[rank_data$group == "b"] <- NA_real_
rank_fit <- crc_fit(rank_data, ~ source_1 + source_2 + source_3,
  ~ source_1 + source_2 + source_3 + group,
  outcome = "outcome", outcome_formula = ~ group)
expect_error(misclassCRC:::maximize_outcome_crc(
  rank_fit$model_matrices, rank_fit$initialization$state_weights, "ztnegbin"),
  pattern = "rank deficient")
