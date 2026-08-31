library(misclassCRC)

data <- data.frame(
  source_1 = rep(c(1L, 0L, 0L, 1L, 1L, 0L), 4L),
  source_2 = rep(c(0L, 1L, 0L, 1L, 0L, 1L), 4L),
  source_3 = rep(c(0L, 0L, 1L, 0L, 1L, 1L), 4L),
  group = factor(rep(c("a", "b"), 12L)),
  lower = rep(c(1, NA, 2, 3, 1, 4), 4L),
  upper = rep(c(1, NA, Inf, 5, 1, 4), 4L)
)
fits <- lapply(c("ztpois", "ztnegbin"), function(distribution) {
  set.seed(42)
  suppressWarnings(crc_fit(
    data, ~ source_1 + source_2 + source_3,
    ~ source_1 + source_2 + source_3 + group + .latent + source_1:.latent,
    outcome = c(lower = "lower", upper = "upper"),
    outcome_formula = ~ group + .latent, outcome_dist = distribution, latent_classes = 2L
  ))
})
for (fit in fits) {
  prediction <- predict(fit)
  grouped <- predict(fit, ~ group + .latent)
  expect_true(is.list(prediction))
  expect_inherits(prediction, "crc_pred")
  expect_equal(names(prediction), c("population_size", "outcome_total", "by_groups"))
  expect_true(all(is.finite(c(prediction$population_size, prediction$outcome_total))))
  expect_equal(prediction$population_size, sum(fit$em_fit$state_weights) +
    sum(fit$em_fit$capture_fit$fitted_mean[fit$model_matrices$capture$unobserved]))
  expect_equal(sum(grouped$by_groups$population_size), prediction$population_size)
  expect_equal(sum(grouped$by_groups$outcome_total), prediction$outcome_total)
  expect_equal(nrow(grouped$by_groups), 4L)
}
expect_equal(sum(predict(fits[[1L]], ~ group)$by_groups$population_size),
  predict(fits[[1L]])$population_size)
expect_equal(sum(predict(fits[[1L]], ~ .latent)$by_groups$outcome_total),
  predict(fits[[1L]])$outcome_total)
limited <- suppressWarnings(crc_fit(data, ~ source_1 + source_2 + source_3,
  ~ source_1 + source_2 + source_3 + group, control = list(em = list(max_iter = 1L))))
expect_false(limited$em_fit$converged)
expect_true(is.finite(predict(limited)$population_size))
expect_null(predict(limited)$outcome_total)
expect_false("outcome_total" %in% names(predict(limited, ~ group)$by_groups))
edge <- fits[[1L]]
edge$model_matrices$states$outcome[1L] <- list(1, Inf, "right_censored")
edge$model_matrices$states$outcome[2L] <- list(500, 501, "interval_censored")
edge$em_fit$outcome_fit$fitted_mean[1:2] <- 10
expect_true(is.finite(predict(edge)$outcome_total))
expect_error(predict(fits[[1L]], ~ group:.latent), pattern = "untransformed")
expect_error(predict(fits[[1L]], ~ unknown), pattern = "Invalid grouping")
expect_error(predict(fits[[1L]], ~ 1), pattern = "at least one")
expect_error(predict(structure(list(fitted = FALSE), class = "crc_fit")), pattern = "fitted `crc_fit`")
