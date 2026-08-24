library(misclassCRC)

data <- data.frame(
  source_1 = c(1L, 0L, 0L, 1L),
  source_2 = c(0L, 1L, 1L, 0L),
  source_3 = c(0L, 0L, 1L, 1L)
)

fit <- crc_fit(
  data = data,
  captures = ~ source_1 + source_2 + source_3,
  capture_formula = ~ source_1 + source_2 + source_3
)
matrices <- fit$model_matrices

expect_inherits(fit, "crcfit")
expect_inherits(matrices, "crc_model_matrices")
expect_equal(dim(matrices$capture$matrix), c(8L, 4L))
expect_equal(sum(matrices$capture$unobserved), 1L)
expect_equal(nrow(matrices$states$data), nrow(data))
expect_equal(matrices$states$misclass_probability, rep(1, nrow(data)))
expect_null(matrices$outcome_matrix)
expect_false(anyNA(matrices$capture_index))

key <- c("source_1", "source_2", "source_3", ".latent")
indexed_cells <- matrices$capture$cells[
  matrices$capture_index, key, with = FALSE
]
expect_equal(indexed_cells, matrices$states$data[, key, with = FALSE])

data$group <- factor(c("b", "a", "b", "a"), levels = c("b", "a"))
data$outcome <- c(2, NA, 4, 5)
fit <- crc_fit(
  data = data,
  captures = ~ source_1 + source_2 + source_3,
  capture_formula = ~ source_1 + source_2 + source_3 + group + .latent,
  outcome = "outcome",
  outcome_formula = ~ group + .latent,
  latent_classes = 2L
)
matrices <- fit$model_matrices

expect_equal(dim(matrices$capture$matrix), c(32L, 6L))
expect_equal(nrow(matrices$states$data), 8L)
expect_equal(dim(matrices$outcome_matrix), c(8L, 3L))
expect_equal(levels(matrices$capture$cells$group), c("b", "a"))
expect_equal(levels(matrices$states$data$group), c("b", "a"))
expect_equal(matrices$states$outcome$lower, data$outcome[matrices$states$observation_id])
expect_equal(length(matrices$capture_index), nrow(matrices$states$data))
