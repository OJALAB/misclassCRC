library(misclassCRC)

data <- data.frame(
  source_1 = c(1L, 0L, 0L, 1L),
  source_2 = c(0L, 1L, 1L, 0L),
  source_3 = c(0L, 0L, 1L, 1L),
  group = c("B", "A", "B", "A"),
  contract = c("perm", "temp", "perm", "temp")
)
group_matrix <- matrix(
  c(0.8, 0.2, 0.1, 0.9), nrow = 2L, byrow = TRUE,
  dimnames = list(c("A", "B"), c("A", "B"))
)
contract_matrix <- matrix(
  c(0.7, 0.3, 0.25, 0.75), nrow = 2L, byrow = TRUE,
  dimnames = list(c("temp", "perm"), c("temp", "perm"))
)

fit <- crc_fit(
  data = data,
  captures = ~ source_1 + source_2 + source_3,
  capture_formula = ~ source_1 + source_2 + source_3 +
    group + contract + .latent,
  misclass = list(
    group = list(matrix = group_matrix, true_if = ~ source_1),
    contract = list(matrix = contract_matrix, true_if = ~ source_3)
  ),
  latent_classes = 2L
)
matrices <- fit$model_matrices
states <- matrices$states

expect_equal(dim(matrices$capture$matrix), c(64L, 7L))
expect_equal(nrow(states$data), 32L)
expect_equal(levels(states$data$group), c("A", "B"))
expect_equal(levels(states$data$contract), c("temp", "perm"))
expect_equal(nrow(states$observed_misclass), nrow(states$data))
expect_equal(length(states$observation_id), nrow(states$data))
expect_equal(length(states$misclass_probability), nrow(states$data))

select_state <- function(id, group, contract) {
  states$observation_id == id & states$data$group == group &
    states$data$contract == contract
}
expect_equal(unique(states$misclass_probability[select_state(2L, "A", "temp")]), 0.56)
expect_equal(unique(states$misclass_probability[select_state(2L, "B", "perm")]), 0.025)
expect_true(all(states$misclass_probability[
  states$observation_id == 1L & states$data$group == "A"
] == 0))
expect_true(all(states$misclass_probability[
  states$observation_id == 3L & states$data$contract == "temp"
] == 0))
expect_null(matrices$outcome_matrix)
expect_false(anyNA(matrices$capture_index))
