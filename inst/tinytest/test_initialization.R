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
prepare_call <- function() misclassCRC:::prepare_crc(
  data = data,
  captures = ~ source_1 + source_2 + source_3,
  capture_formula = ~ source_1 + source_2 + source_3 + group + .latent,
  outcome = NULL,
  outcome_formula = NULL,
  outcome_dist = "ztnegbin",
  misclass = list(group = list(matrix = group_matrix, true_if = ~ source_1)),
  latent_classes = 2L,
  control = list(init_alpha = 12),
  verbose = FALSE
)

set.seed(42)
preparation <- prepare_call()
initialization <- preparation$initialization
matrices <- preparation$model_matrices
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
expect_equal(prepare_call()$initialization, initialization)

preparation_one <- misclassCRC:::prepare_crc(
  data, ~ source_1 + source_2 + source_3,
  ~ source_1 + source_2 + source_3,
  outcome = NULL, outcome_formula = NULL, outcome_dist = "ztnegbin",
  misclass = NULL, latent_classes = 1L, control = NULL, verbose = FALSE
)
expect_equal(preparation_one$initialization$state_weights, rep(1, nrow(data)))
expect_equal(preparation_one$initialization$init_alpha, 20)
expect_error(misclassCRC:::prepare_crc(
  data, ~ source_1 + source_2 + source_3,
  ~ source_1 + source_2 + source_3,
  outcome = NULL, outcome_formula = NULL, outcome_dist = "ztnegbin",
  misclass = NULL, latent_classes = 1L,
  control = list(init_alpha = 0), verbose = FALSE
),
  pattern = "finite positive number")
