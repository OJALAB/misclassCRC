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
fit_call <- function() crc_fit(
  data = data,
  captures = ~ source_1 + source_2 + source_3,
  capture_formula = ~ source_1 + source_2 + source_3 + group + .latent,
  misclass = list(group = list(matrix = group_matrix, true_if = ~ source_1)),
  latent_classes = 2L,
  control = list(init_alpha = 12)
)

set.seed(42)
fit <- fit_call()
initialization <- fit$initialization
matrices <- fit$model_matrices
state_sums <- tapply(initialization$state_weights,
  matrices$states$observation_id, sum)
expected_counts <- numeric(nrow(matrices$capture$cells))
for (i in seq_along(initialization$state_weights)) {
  cell <- matrices$capture_index[[i]]
  expected_counts[[cell]] <- expected_counts[[cell]] + initialization$state_weights[[i]]
}
expected_counts[matrices$capture$unobserved] <- 1
expect_inherits(initialization, "crc_initialization")
expect_equal(length(initialization$state_weights), nrow(matrices$states$data))
expect_true(all(is.finite(initialization$state_weights)))
expect_true(all(initialization$state_weights >= 0))
expect_equal(as.numeric(state_sums), rep(1, nrow(data)))
expect_equal(initialization$cell_counts, expected_counts)
expect_equal(initialization$init_alpha, 12)
set.seed(42)
expect_equal(fit_call()$initialization, initialization)

fit_one <- crc_fit(data, ~ source_1 + source_2 + source_3,
  ~ source_1 + source_2 + source_3)
expect_equal(fit_one$initialization$state_weights, rep(1, nrow(data)))
expect_equal(fit_one$initialization$init_alpha, 20)
expect_error(crc_fit(data, ~ source_1 + source_2 + source_3,
  ~ source_1 + source_2 + source_3, control = list(init_alpha = 0)),
  pattern = "finite positive number")
