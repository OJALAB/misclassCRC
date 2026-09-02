library(data.table)
library(misclassCRC)
library(futurize)
library(progressify)

source("simulations/functions.R")
load("simulations/data.RData")

population <- as.data.table(data)[
  substr(poz_kodZawodu, 1L, 1L) != "0",
  .(
    code = factor(substr(poz_kodZawodu, 1L, 1L), levels = as.character(1:9)),
    vac = as.numeric(poz_lWolnychMiejsc)
  )
]
population[, code := relevel(code, ref = "2")]
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

n_rep <- 100L
workers <- 10L

set.seed(123L)

plan(multisession, workers = workers)

results <- rbindlist(lapply(1:n_rep, function(i) {

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

  data.table(
    population_size = overall_prediction$population_size,
    outcome_total = overall_prediction$outcome_total
  )

}) |> progressify() |> futurize(seed = TRUE))

results[, `:=` (
  true_population_size = rep(nrow(population), n_rep),
  true_outcome_total = rep(sum(population[["vac"]]), n_rep)
)]

estimands <- data.table(
  estimand = c("population_size", "outcome_total"),
  estimate = c("population_size", "outcome_total"),
  truth = c("true_population_size", "true_outcome_total")
)

evaluation <- rbindlist(lapply(seq_len(nrow(estimands)), function(i) {
  estimate <- results[[estimands$estimate[i]]]
  truth <- results[[estimands$truth[i]]]

  error <- estimate - truth

  data.table(
    estimand = estimands$estimand[i],
    mean_estimate = mean(estimate),
    MB = mean(error),
    MRB = 100 * mean(error / truth),
    RMSE = sqrt(mean(error^2))
  )
}))
