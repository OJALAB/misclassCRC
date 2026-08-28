library(misclassCRC)

data <- data.frame(
  source_1 = rep(c(1L, 0L, 0L, 1L, 1L, 0L), 4L),
  source_2 = rep(c(0L, 1L, 0L, 1L, 0L, 1L), 4L),
  source_3 = rep(c(0L, 0L, 1L, 0L, 1L, 1L), 4L),
  group = factor(rep(c("a", "b"), 12L)),
  outcome = rep(c(1, 2, 3, 10, 2, 14), 4L)
)

fit <- crc_fit(
  data,
  captures = ~ source_1 + source_2 + source_3,
  capture_formula = ~ source_1 + source_2 + source_3 + group,
  outcome = "outcome",
  outcome_formula = ~ group,
  outcome_dist = "ztpois"
)

expect_inherits(fit, "crc_fit")
expect_true(fit$fitted)
expect_inherits(fit$em_fit, "crc_em_fit")
expect_true(fit$em_fit$converged)
expect_inherits(fit$em_fit$capture_fit, "crc_capture_mstep")
expect_inherits(fit$em_fit$outcome_fit, "crc_outcome_mstep")
expect_equal(fit$outcome_dist, "ztpois")
expect_equal(fit$call[[1L]], quote(crc_fit))

without_outcome <- crc_fit(
  data,
  captures = ~ source_1 + source_2 + source_3,
  capture_formula = ~ source_1 + source_2 + source_3 + group
)

expect_true(without_outcome$em_fit$converged)
expect_null(without_outcome$em_fit$outcome_fit)
expect_null(without_outcome$outcome_model)

expect_warning(
  limited <- crc_fit(
    data,
    captures = ~ source_1 + source_2 + source_3,
    capture_formula = ~ source_1 + source_2 + source_3 + group,
    control = list(em = list(max_iter = 1L))
  ),
  pattern = "did not converge after 1 iteration"
)
expect_false(limited$em_fit$converged)
expect_equal(limited$em_fit$iterations, 1L)
