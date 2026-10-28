library(misclassCRC)

data <- data.frame(
  source_1 = rep(c(1L, 0L, 0L, 1L, 1L, 0L), 4L),
  source_2 = rep(c(0L, 1L, 0L, 1L, 0L, 1L), 4L),
  source_3 = rep(c(0L, 0L, 1L, 0L, 1L, 1L), 4L),
  group = factor(rep(c("a", "b"), 12L)),
  outcome = rep(c(1, 2, 3, 10, 2, 14), 4L),
  lower = rep(c(1, NA, 2, 3, 1, 4), 4L),
  upper = rep(c(1, NA, Inf, 5, 1, 4), 4L)
)
captures <- ~ source_1 + source_2 + source_3

has_line <- function(output, pattern) any(grepl(pattern, output))

set.seed(42)
nb_fit <- suppressWarnings(crc_fit(
  data, captures,
  ~ source_1 + source_2 + source_3 + group + .latent + source_1:.latent,
  outcome = c(lower = "lower", upper = "upper"),
  outcome_formula = ~ group + .latent, outcome_dist = "ztnegbin",
  latent_classes = 2L
))
pois_fit <- crc_fit(
  data, captures, ~ source_1 + source_2 + source_3 + group,
  outcome = "outcome", outcome_formula = ~ group, outcome_dist = "ztpois"
)
without_outcome <- crc_fit(data, captures, ~ source_1 + source_2 + source_3 + group)
limited <- suppressWarnings(crc_fit(
  data, captures, ~ source_1 + source_2 + source_3 + group,
  control = list(em = list(max_iter = 1L))
))
group_matrix <- matrix(
  c(0.9, 0.1, 0.2, 0.8), nrow = 2L, byrow = TRUE,
  dimnames = list(c("a", "b"), c("a", "b"))
)
misclass_fit <- suppressWarnings(crc_fit(
  data, captures, ~ source_1 + source_2 + source_3 + group,
  misclass = list(group = list(matrix = group_matrix, true_if = ~ source_1))
))

# print.crc_fit()

printed <- capture.output(visible <- withVisible(print(nb_fit)))
expect_false(visible$visible)
expect_identical(visible$value, nb_fit)
expect_true(has_line(printed, "^Capture-recapture model fit$"))
expect_true(has_line(printed, "^Observed units: +24$"))
expect_true(has_line(printed, "^Capture sources: +source_1, source_2, source_3$"))
expect_true(has_line(printed, "^Latent classes: +2$"))
expect_true(has_line(printed, "^Misclassified variables: +none$"))
expect_true(has_line(printed, "^Outcome distribution: +zero-truncated negative binomial$"))
expect_true(has_line(capture.output(print(pois_fit)), "zero-truncated Poisson$"))
expect_true(has_line(capture.output(print(pois_fit)), "^EM algorithm: +converged after"))
expect_true(has_line(capture.output(print(without_outcome)), "^Outcome distribution: +none$"))
expect_true(has_line(
  capture.output(print(limited)), "^EM algorithm: +not converged after 1 iteration$"
))
expect_true(has_line(capture.output(print(misclass_fit)), "^Misclassified variables: +group$"))

# summary.crc_fit()

nb_summary <- summary(nb_fit)
expect_inherits(nb_summary, "summary.crc_fit")
expect_equal(
  names(nb_summary),
  c("call", "model", "misclassification", "capture", "outcome",
    "convergence", "converged", "iterations")
)
expect_equal(nb_summary$call, nb_fit$call)
expect_equal(nb_summary$model[["Latent classes"]], 2L)
expect_equal(
  nb_summary$model[["Capture cells"]], nrow(nb_fit$model_matrices$capture$cells)
)
expect_equal(nb_summary$capture$coefficients, nb_fit$em_fit$capture_fit$coefficients)
expect_equal(nb_summary$outcome$coefficients, nb_fit$em_fit$outcome_fit$coefficients)
expect_true(is.finite(nb_summary$outcome$sigma) && nb_summary$outcome$sigma > 0)
expect_null(summary(pois_fit)$outcome$sigma)
expect_null(summary(without_outcome)$outcome)
expect_equal(length(nb_summary$misclassification), 0L)
expect_equal(
  summary(misclass_fit)$misclassification,
  c(group = "2 levels; true for source_1")
)
expect_false(summary(limited)$converged)
expect_equal(summary(limited)$iterations, 1L)

# print.summary.crc_fit()

printed <- capture.output(visible <- withVisible(print(nb_summary)))
expect_false(visible$visible)
expect_identical(visible$value, nb_summary)
expect_true(has_line(printed, "^Capture-recapture model summary$"))
expect_true(has_line(printed, "^Capture model:$"))
expect_true(has_line(printed, "^Outcome model:$"))
expect_true(has_line(printed, "Dispersion \\(sigma\\):"))
expect_false(has_line(printed, "^Misclassification:$"))
pois_printed <- capture.output(print(summary(pois_fit)))
expect_true(has_line(pois_printed, "^Outcome model:$"))
expect_false(has_line(pois_printed, "Dispersion"))
expect_false(has_line(capture.output(print(summary(without_outcome))), "^Outcome model:$"))
expect_true(has_line(
  capture.output(print(summary(misclass_fit))), "^Misclassification:$"
))
expect_true(has_line(capture.output(print(summary(limited))), "Status: +not converged$"))
expect_false(identical(
  capture.output(print(nb_summary, digits = 3L)),
  capture.output(print(nb_summary, digits = 7L))
))

# print.crc_pred()

prediction <- predict(nb_fit)
grouped <- predict(nb_fit, ~ group + .latent)
no_outcome <- predict(without_outcome, ~ group)

printed <- capture.output(visible <- withVisible(print(prediction)))
expect_false(visible$visible)
expect_identical(visible$value, prediction)
expect_true(has_line(printed, "^Capture-recapture prediction$"))
expect_true(has_line(printed, "^Population size:"))
expect_true(has_line(printed, "^Outcome total:"))
expect_false(has_line(printed, "Grouping variables"))
expect_false(has_line(printed, "not modeled"))
printed <- capture.output(print(grouped))
expect_true(has_line(printed, "^Grouping variables: +group, \\.latent$"))
expect_true(has_line(printed, "^Groups: +4$"))
expect_true(has_line(printed, "^Grouped estimates:$"))
expect_true(has_line(capture.output(print(no_outcome)), "^Outcome total: +not modeled$"))
large <- prediction
large$population_size <- 1234567.891
large$outcome_total <- 9876543.21
printed <- capture.output(print(large))
expect_true(has_line(printed, "^Population size: +1,234,568$"))
expect_true(has_line(printed, "^Outcome total: +9,876,543$"))
expect_false(has_line(printed, "e\\+"))

# summary.crc_pred()

prediction_summary <- summary(prediction)
expect_inherits(prediction_summary, "summary.crc_pred")
expect_equal(names(prediction_summary), c("overall", "groups", "by_groups"))
expect_equal(
  names(prediction_summary$overall),
  c("population_size", "outcome_total", "mean_outcome")
)
expect_equal(prediction_summary$overall[["population_size"]], prediction$population_size)
expect_equal(
  prediction_summary$overall[["mean_outcome"]],
  prediction$outcome_total / prediction$population_size
)
expect_identical(prediction_summary$groups, character())
expect_null(prediction_summary$by_groups)

grouped_summary <- summary(grouped)
by_groups <- grouped_summary$by_groups
expect_identical(grouped_summary$groups, c("group", ".latent"))
expect_equal(
  names(by_groups),
  c("group", ".latent", "population_size", "population_share",
    "outcome_total", "outcome_share", "mean_outcome")
)
expect_equal(sum(by_groups$population_share), 1)
expect_equal(sum(by_groups$outcome_share), 1)
expect_equal(by_groups$population_size, grouped$by_groups$population_size)
expect_equal(by_groups$mean_outcome, by_groups$outcome_total / by_groups$population_size)
expect_equal(
  names(grouped$by_groups),
  c("group", ".latent", "population_size", "outcome_total")
)

no_outcome_summary <- summary(no_outcome)
expect_equal(names(no_outcome_summary$overall), "population_size")
expect_equal(
  names(no_outcome_summary$by_groups),
  c("group", "population_size", "population_share")
)

# print.summary.crc_pred()

printed <- capture.output(visible <- withVisible(print(grouped_summary)))
expect_false(visible$visible)
expect_identical(visible$value, grouped_summary)
expect_true(has_line(printed, "^Capture-recapture prediction summary$"))
expect_true(has_line(printed, "^Overall estimates:$"))
expect_true(has_line(printed, "Mean outcome per unit:"))
expect_true(has_line(printed, "^Grouped estimates:$"))
expect_true(has_line(printed, "Grouping variables: +group, \\.latent$"))
expect_true(has_line(printed, "Groups: +4$"))
expect_false(has_line(capture.output(print(prediction_summary)), "^Grouped estimates:$"))
printed <- capture.output(print(no_outcome_summary))
expect_true(has_line(printed, "Outcome total: +not modeled$"))
expect_false(has_line(printed, "Mean outcome per unit"))
expect_false(identical(
  capture.output(print(grouped_summary, digits = 3L)),
  capture.output(print(grouped_summary, digits = 7L))
))
