library(misclassCRC)

data <- data.frame(
  source_1 = c(1L, 0L, 0L, 1L),
  source_2 = c(0L, 1L, 1L, 0L),
  source_3 = c(0L, 0L, 1L, 1L),
  group = c("B", "A", "B", "A")
)
group_matrix <- matrix(
  c(0.8, 0.2, 0.1, 0.9), nrow = 2L, byrow = TRUE,
  dimnames = list(c("A", "B"), c("A", "B"))
)
fit <- crc_fit(
  data = data,
  captures = ~ source_1 + source_2 + source_3,
  capture_formula = ~ source_1 + source_2 + source_3 + group + .latent,
  misclass = list(group = list(matrix = group_matrix, true_if = ~ source_1)),
  latent_classes = 2L
)
matrices <- fit$model_matrices
n_states <- nrow(matrices$states$data)
n_cells <- nrow(matrices$capture$cells)
capture_mean <- exp(seq(-1, 1, length.out = n_cells))
outcome_loglik <- rep(c(-0.25, -0.75), length.out = n_states)

result <- misclassCRC:::expectation_crc(matrices, capture_mean, outcome_loglik)
raw <- capture_mean[matrices$capture_index] *
  matrices$states$misclass_probability * exp(outcome_loglik)
expected <- raw / ave(raw, matrices$states$observation_id, FUN = sum)

expect_inherits(result, "crc_estep")
expect_equal(result$state_weights, expected)
expect_equal(length(result$state_weights), n_states)
expect_equal(length(result$cell_counts), n_cells)
expect_equal(as.numeric(tapply(result$state_weights,
  matrices$states$observation_id, sum)), rep(1, nrow(data)))
expect_true(all(result$state_weights[matrices$states$misclass_probability == 0] == 0))
expect_equal(sum(result$cell_counts[!matrices$capture$unobserved]), nrow(data))
expect_equal(result$cell_counts[matrices$capture$unobserved],
  capture_mean[matrices$capture$unobserved])
expect_equal(
  misclassCRC:::expectation_crc(matrices, capture_mean)$state_weights,
  misclassCRC:::expectation_crc(
    matrices,
    capture_mean,
    numeric(n_states)
  )$state_weights
)

stable <- misclassCRC:::expectation_crc(
  matrices,
  capture_mean,
  outcome_loglik + 700
)
expect_true(all(is.finite(stable$state_weights)))
expect_equal(stable$state_weights, result$state_weights)
expect_error(
  misclassCRC:::expectation_crc(matrices, capture_mean[-1]),
  pattern = "not aligned"
)
expect_error(
  misclassCRC:::expectation_crc(matrices, replace(capture_mean, 1L, 0)),
  pattern = "non-positive"
)
expect_error(
  misclassCRC:::expectation_crc(matrices, capture_mean, outcome_loglik[-1]),
  pattern = "invalid"
)
expect_error(
  misclassCRC:::expectation_crc(matrices, capture_mean, rep(-Inf, n_states)),
  pattern = "no admissible"
)
