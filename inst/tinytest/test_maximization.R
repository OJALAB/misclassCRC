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
