generate_iaos_data <- function(
  population, confusion_matrix, p_easy, p_hard,
  prob_hard, validation_size, missing_fraction,
  censoring_fraction, censoring_threshold
) {
  population <- data.table::copy(data.table::as.data.table(population))
  if (!all(c("code", "vac") %in% names(population))) {
    stop("`population` must contain `code` and `vac` columns.", call. = FALSE)
  }
  code_levels <- c("2", "1", "3", "4", "5", "6", "7", "8", "9")
  code_values <- as.character(population[["code"]])
  if (anyNA(code_values) || any(!code_values %in% as.character(1:9))) {
    stop("`code` must contain only categories 1 through 9.", call. = FALSE)
  }
  population[, code := factor(code_values, levels = code_levels)]
  if (!is.numeric(population[["vac"]]) || anyNA(population[["vac"]]) ||
      any(!is.finite(population[["vac"]])) || any(population[["vac"]] < 1) ||
      any(population[["vac"]] != floor(population[["vac"]]))) {
    stop("`vac` must be positive-integer numeric.", call. = FALSE)
  }
  population[, vac := as.numeric(vac)]

  if (!is.matrix(confusion_matrix) || !is.numeric(confusion_matrix) ||
      !identical(dim(confusion_matrix), c(9L, 9L)) ||
      is.null(rownames(confusion_matrix)) || is.null(colnames(confusion_matrix)) ||
      !setequal(rownames(confusion_matrix), as.character(1:9)) ||
      !setequal(colnames(confusion_matrix), as.character(1:9))) {
    stop("`confusion_matrix` must be a named numeric 9 by 9 matrix.", call. = FALSE)
  }
  confusion_matrix <- confusion_matrix[as.character(1:9), as.character(1:9)]
  row_totals <- rowSums(confusion_matrix)
  if (anyNA(confusion_matrix) || any(!is.finite(confusion_matrix)) ||
      any(confusion_matrix < 0) || any(!is.finite(row_totals)) ||
      any(abs(row_totals - 1) > sqrt(.Machine$double.eps))) {
    stop("Each confusion-matrix row must contain probabilities summing to one.", call. = FALSE)
  }
  check_probabilities <- function(x, n, label) {
    if (!is.numeric(x) || length(x) != n || anyNA(x) ||
        any(!is.finite(x)) || any(x < 0 | x > 1)) {
      stop(sprintf("`%s` must contain %d valid probabilities.", label, n), call. = FALSE)
    }
  }
  check_probabilities(p_easy, 3L, "p_easy")
  check_probabilities(p_hard, 3L, "p_hard")
  check_probabilities(prob_hard, 9L, "prob_hard")
  check_probabilities(missing_fraction, 1L, "missing_fraction")
  check_probabilities(censoring_fraction, 1L, "censoring_fraction")
  if (missing_fraction + censoring_fraction > 1) {
    stop("Missing and censoring fractions must sum to at most one.", call. = FALSE)
  }
  if (length(validation_size) != 1L || is.na(validation_size) ||
      validation_size < 1 || validation_size != as.integer(validation_size)) {
    stop("`validation_size` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(censoring_threshold) || length(censoring_threshold) != 1L ||
      is.na(censoring_threshold) || !is.finite(censoring_threshold) ||
      censoring_threshold < 1) {
    stop("`censoring_threshold` must be a finite positive number.", call. = FALSE)
  }

  generate_labels <- function(x) {
    x[, code_ml := {
      group <- as.character(code[1L])
      as.character(sample(
        colnames(confusion_matrix), .N, replace = TRUE,
        prob = confusion_matrix[group, ]
      ))
    }, by = .(code)]
    x[, code_ml := factor(code_ml, levels = code_levels)]
    x
  }
  population <- generate_labels(population)
  code_distribution <- prop.table(table(population[["code"]]))
  validation_data <- data.table::data.table(code = factor(
    sample(names(code_distribution), validation_size, replace = TRUE,
      prob = as.numeric(code_distribution)),
    levels = code_levels
  ))
  validation_data <- generate_labels(validation_data)
  validation_counts <- table(validation_data[["code"]], validation_data[["code_ml"]])
  if (any(rowSums(validation_counts) == 0L)) {
    stop("The validation sample contains an empty true-code row.", call. = FALSE)
  }
  estimated_matrix <- unclass(prop.table(validation_counts, margin = 1L))
  storage.mode(estimated_matrix) <- "double"

  hard_parameters <- data.table::data.table(
    code = factor(as.character(1:9), levels = code_levels),
    probability_hard = prob_hard
  )
  population <- merge(population, hard_parameters, by = "code", sort = FALSE)
  population[, is_hard := rbinom(.N, 1L, probability_hard)]
  population[, `:=`(
    probability_1 = ifelse(is_hard == 1L, p_hard[1L], p_easy[1L]),
    probability_2 = ifelse(is_hard == 1L, p_hard[2L], p_easy[2L]),
    probability_3 = ifelse(is_hard == 1L, p_hard[3L], p_easy[3L])
  )]
  population[, I1 := rbinom(.N, 1L, probability_1)]
  population[, probability_2 := pmax(probability_2 - 0.05 * I1, 0)]
  population[, I2 := rbinom(.N, 1L, probability_2)]
  population[, probability_3 := pmax(probability_3 - 0.05 * I1, 0)]
  population[, I3 := rbinom(.N, 1L, probability_3)]

  observed <- data.table::copy(population[I1 + I2 + I3 > 0L])
  observed[, observed_code := code_ml]
  observed[I1 == 1L, observed_code := code]
  observed[, `:=`(vac_lower = vac, vac_upper = vac)]
  observed[I1 == 0L, bound_draw := runif(.N)]
  observed[I1 == 0L & bound_draw <= missing_fraction,
    `:=`(vac_lower = NA_real_, vac_upper = NA_real_)]
  observed[I1 == 0L & !is.na(vac_lower) & vac >= censoring_threshold &
      bound_draw <= missing_fraction + censoring_fraction,
    `:=`(vac_lower = as.numeric(censoring_threshold), vac_upper = Inf)]
  if (any(xor(is.na(observed[["vac_lower"]]), is.na(observed[["vac_upper"]]))) ||
      any(observed[!is.na(vac_lower), vac_lower > vac_upper])) {
    stop("Generated outcome bounds are invalid.", call. = FALSE)
  }
  observed[, code := factor(observed_code, levels = code_levels)]
  observed_data <- observed[, .(
    I1 = as.integer(I1), I2 = as.integer(I2), I3 = as.integer(I3),
    code, vac_lower = as.numeric(vac_lower), vac_upper = as.numeric(vac_upper)
  )]

  truth_overall <- population[, .(
    population_size = .N, outcome_total = sum(vac)
  )]
  truth_by_code <- population[, .(
    population_size = .N, outcome_total = sum(vac)
  ), by = .(code)]
  truth_by_code[, code_order := as.integer(as.character(code))]
  data.table::setorder(truth_by_code, code_order)
  truth_by_code[, code_order := NULL]

  list(
    data = observed_data,
    confusion_matrix = estimated_matrix,
    validation_data = validation_data[, .(code, code_ml)],
    truth = list(
      population = population,
      overall = truth_overall,
      by_code = truth_by_code
    )
  )
}
