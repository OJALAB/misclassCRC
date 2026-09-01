library(data.table)
library(misclassCRC)

source("simulations/functions.R")
load("simulations/data.RData")

population <- data.table::as.data.table(data)[
  substr(poz_kodZawodu, 1L, 1L) != "0",
  .(
    code = factor(substr(poz_kodZawodu, 1L, 1L), levels = as.character(1:9)),
    vac = as.numeric(poz_lWolnychMiejsc)
  )
]
population[, code := relevel(code, ref = "2")]
if (nrow(population) != 10524L || sum(population[["vac"]]) != 21532) {
  stop("The IAOS population must have size 10,524 and vacancy total 21,532.",
    call. = FALSE)
}

confusion_percentages <- c(
  82.8,11.2,1.7,1.7,2.2,0,0.4,0,0, 3.5,88.2,4.8,1.6,1.6,0,0.3,0,0,
  2,13.9,68.8,2.6,4.4,0,6,1.2,1.2, 1.1,10.2,6,71.3,3,0,0,1.1,7.2,
  1.5,4.3,6.3,1,85,0,0.8,0,1.3, 0,14.3,0,0,0,71.4,0,0,14.3,
  0.6,1.4,1.7,0.2,0.4,0,90.7,2.7,2.3, 0,0,1.9,1.2,0.4,0,18.2,77.9,0.4,
  0.5,0,2.7,5,2.7,2.3,9.1,0.9,76.8
)
confusion_matrix <- matrix(confusion_percentages / 100, 9L, 9L, byrow = TRUE,
  dimnames = list(as.character(1:9), as.character(1:9)))
confusion_matrix <- confusion_matrix / rowSums(confusion_matrix)

iaos_parameters <- list(
  seed = 321L,
  p_easy = c(0.8, 0.2, 0.1),
  p_hard = c(0.1, 0.8, 0.6),
  prob_hard = c(0.4, 0.3, 0.2, 0.3, 0.4, 0.3, 0.2, 0.3, 0.4),
  validation_size = 1000L,
  missing_fraction = 0.25,
  censoring_fraction = 0.25,
  censoring_threshold = 2,
  em_max_iter = 2000L,
  em_tolerance = 5e-4
)

set.seed(iaos_parameters$seed)
generated <- generate_iaos_data(
  population = population,
  confusion_matrix = confusion_matrix,
  p_easy = iaos_parameters$p_easy,
  p_hard = iaos_parameters$p_hard,
  prob_hard = iaos_parameters$prob_hard,
  validation_size = iaos_parameters$validation_size,
  missing_fraction = iaos_parameters$missing_fraction,
  censoring_fraction = iaos_parameters$censoring_fraction,
  censoring_threshold = iaos_parameters$censoring_threshold
)

fit <- crc_fit(
  data = generated$data,
  captures = ~ I1 + I2 + I3,
  capture_formula = ~ (I1 + I2 + I3) * .latent +
    code * .latent + I1:I2 + I1:I3,
  outcome = c(lower = "vac_lower", upper = "vac_upper"),
  outcome_formula = ~ code,
  outcome_dist = "ztnegbin",
  misclass = list(code = list(
    matrix = generated$confusion_matrix, true_if = ~ I1
  )),
  latent_classes = 2L,
  control = list(em = list(max_iter = iaos_parameters$em_max_iter, tolerance = iaos_parameters$em_tolerance))
)

overall_prediction <- predict(fit)
by_code_prediction <- predict(fit, by = ~ code)
overall_comparison <- data.table(
  measure = c("population_size", "outcome_total"),
  truth = c(generated$truth$overall$population_size,
    generated$truth$overall$outcome_total),
  estimate = c(overall_prediction$population_size,
    overall_prediction$outcome_total)
)
overall_comparison[, difference := estimate - truth]

by_code_estimates <- copy(by_code_prediction$by_groups)
by_code_comparison <- merge(
  generated$truth$by_code,
  by_code_estimates,
  by = "code",
  suffixes = c("_truth", "_estimate"),
  all = TRUE,
  sort = FALSE
)
by_code_comparison[, `:=`(
  population_size_difference = population_size_estimate - population_size_truth,
  outcome_total_difference = outcome_total_estimate - outcome_total_truth,
  code_order = as.integer(as.character(code))
)]
setorder(by_code_comparison, code_order)
by_code_comparison[, code_order := NULL]

bound_status_counts <- generated$data[, .(
  exact = sum(!is.na(vac_lower) & is.finite(vac_upper) & vac_lower == vac_upper),
  missing = sum(is.na(vac_lower) & is.na(vac_upper)),
  right_censored = sum(!is.na(vac_lower) & is.infinite(vac_upper))
)]
generated_summary <- list(
  full_population_size = nrow(generated$truth$population),
  observed_population_size = nrow(generated$data),
  all_zero_count = sum(generated$truth$population$I1 +
    generated$truth$population$I2 + generated$truth$population$I3 == 0L),
  capture_totals = colSums(generated$data[, .(I1, I2, I3)]),
  bound_status_counts = bound_status_counts
)
convergence <- list(
  converged = fit$em_fit$converged,
  iterations = fit$em_fit$iterations,
  maximum = fit$em_fit$convergence$maximum,
  changes = fit$em_fit$convergence$changes
)

iaos_result <- list(
  seed = iaos_parameters$seed,
  parameters = iaos_parameters,
  generated_summary = generated_summary,
  confusion_matrix = generated$confusion_matrix,
  fit = fit,
  convergence = convergence,
  truth = generated$truth,
  predictions = list(overall = overall_prediction, by_code = by_code_prediction),
  comparisons = list(overall = overall_comparison, by_code = by_code_comparison)
)

print(overall_comparison)
print(by_code_comparison)
