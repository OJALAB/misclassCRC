#' @import data.table
#' @importFrom stats terms
#' 
#' @noRd
parse_captures <- function(captures, data) {

  if (
    !inherits(captures, "formula") ||
      length(captures) != 2L
  ) {
    stop(
      "`captures` must be a one-sided formula.",
      call. = FALSE
    )
  }

  capture_terms_object <- terms(captures)
  capture_terms <- attr(capture_terms_object, "term.labels")
  capture_names <- all.vars(captures)

  if (
    !identical(capture_terms, capture_names) ||
      attr(capture_terms_object, "intercept") != 1L
  ) {
    stop(
      "`captures` must contain only capture-variable names joined with `+`.",
      call. = FALSE
    )
  }

  if (length(capture_names) < 2L) {
    stop(
      "`captures` must identify at least two capture variables.",
      call. = FALSE
    )
  }

  if (length(capture_names) == 2L) {
    warning(
      "With two capture sources, the unobserved cell is identifiable only ",
      "under restrictive assumptions. Source interactions generally cannot ",
      "be estimated freely.",
      call. = FALSE
    )
  }

  missing_names <- setdiff(capture_names, names(data))

  if (length(missing_names) > 0L) {
    stop(
      "The following capture variables are not present in `data`: ",
      paste0("`", missing_names, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  capture_data <- as.data.table(data)[, capture_names, with = FALSE]

  for (var in capture_names) {
    x <- capture_data[[var]]

    if (anyNA(x)) {
      stop(
        "Capture variable `",
        var,
        "` contains missing values.",
        call. = FALSE
      )
    }

    if (!is.numeric(x) && !is.logical(x)) {
      stop(
        "Capture variable `",
        var,
        "` must be numeric, integer, or logical.",
        call. = FALSE
      )
    }

    if (!all(x %in% c(0L, 1L))) {
      stop(
        "Capture variable `",
        var,
        "` must contain only 0 and 1.",
        call. = FALSE
      )
    }

    if (!any(x == 1L)) {
      stop(
        "Capture variable `",
        var,
        "` does not contain any captured units.",
        call. = FALSE
      )
    }

    if (all(x == 1L)) {
      warning(
        "Capture variable `",
        var,
        "` contains only captured units.",
        call. = FALSE
      )
    }
  }

  capture_data[,
    (capture_names) := lapply(.SD, as.integer),

    .SDcols = capture_names
  ]

  observed <- rowSums(capture_data) > 0L

  if (!all(observed)) {
    stop(
      "Each row of `data` must be observed in at least one capture source.",
      call. = FALSE
    )
  }

  structure(
    list(
      names = capture_names,
      data = capture_data,
      n_sources = length(capture_names)
    ),
    class = "crc_captures"
  )

}

#' @import data.table
#' 
#' @noRd
parse_outcome <- function(outcome, data) {

  if (is.null(outcome)) {
    return(NULL)
  }

  if (is.character(outcome) && length(outcome) == 1L) {
    if (!outcome %in% names(data)) {
      stop(
        "Outcome variable `",
        outcome,
        "` is not present in `data`.",
        call. = FALSE
      )
    }

    x <- data[[outcome]]

    if (!is.numeric(x)) {
      stop(
        "The outcome variable must be numeric.",
        call. = FALSE
      )
    }

    observed_x <- x[!is.na(x)]

    if (
      any(!is.finite(observed_x)) ||
        any(observed_x < 1) ||
        any(observed_x != floor(observed_x))
    ) {
      stop(
        "Observed outcome values must be finite positive integers.",
        call. = FALSE
      )
    }

    if (all(is.na(x))) {
      stop(
        "The outcome variable contains no observed values.",
        call. = FALSE
      )
    }

    return(
      data.table(
        lower = x,
        upper = x,
        status = ifelse(is.na(x), "missing", "exact")
      )
    )
  }

  if (
    !is.character(outcome) ||
      length(outcome) != 2L ||
      is.null(names(outcome)) ||
      anyDuplicated(names(outcome)) > 0L ||
      !setequal(names(outcome), c("lower", "upper"))
  ) {
    stop(
      "`outcome` must be `NULL`, a single column name, or a named vector ",
      "`c(lower = ..., upper = ...)`.",
      call. = FALSE
    )
  }

  missing_names <- setdiff(unname(outcome), names(data))

  if (length(missing_names) > 0L) {
    stop(
      "The following outcome variables are not present in `data`: ",
      paste0("`", missing_names, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  lower <- data[[outcome[["lower"]]]]
  upper <- data[[outcome[["upper"]]]]

  if (!is.numeric(lower) || !is.numeric(upper)) {
    stop(
      "Outcome bounds must be numeric.",
      call. = FALSE
    )
  }

  if (any(lower > upper, na.rm = TRUE)) {
    stop(
      "Outcome lower bounds must not exceed upper bounds.",
      call. = FALSE
    )
  }

  if (any(xor(is.na(lower), is.na(upper)))) {
    stop(
      "Lower and upper outcome bounds must either both be missing or both observed.",
      call. = FALSE
    )
  }

  if (any(!is.finite(lower) & !is.na(lower))) {
    stop(
      "Outcome lower bounds must be finite.",
      call. = FALSE
    )
  }

  if (any(upper == -Inf, na.rm = TRUE)) {
    stop(
      "Outcome upper bounds cannot be `-Inf`.",
      call. = FALSE
    )
  }

  observed_lower <- lower[!is.na(lower)]
  finite_upper <- upper[!is.na(upper) & is.finite(upper)]

  if (
    any(observed_lower < 1) ||
      any(observed_lower != floor(observed_lower)) ||
      any(finite_upper < 1) ||
      any(finite_upper != floor(finite_upper))
  ) {
    stop(
      "Observed outcome bounds must be positive integers.",
      call. = FALSE
    )
  }

  if (all(is.na(lower) & is.na(upper))) {
    stop(
      "The outcome contains no observed or censored values.",
      call. = FALSE
    )
  }

  status <- fcase(
    is.na(lower) & is.na(upper) , "missing"        ,
    lower == upper              , "exact"          ,
    is.infinite(upper)          , "right_censored" ,
    default = "interval_censored"
  )

  data.table(
    lower = lower,
    upper = upper,
    status = status
  )
}

#' @importFrom stats terms
#'
#' @noRd
parse_capture_formula <- function(
  capture_formula,
  data,
  capture_names,
  latent_classes
) {

  if (
    !inherits(capture_formula, "formula") ||
      length(capture_formula) != 2L
  ) {
    stop(
      "`capture_formula` must be a one-sided formula.",
      call. = FALSE
    )
  }

  if (".latent" %in% names(data)) {
    stop(
      "`.latent` is a reserved name and cannot be a column in `data`.",
      call. = FALSE
    )
  }

  terms_object <- terms(capture_formula)
  term_labels <- attr(terms_object, "term.labels")

  if (attr(terms_object, "intercept") != 1L) {
    stop(
      "`capture_formula` must include an intercept.",
      call. = FALSE
    )
  }

  formula_variables <- all.vars(capture_formula)
  ordinary_variables <- setdiff(formula_variables, ".latent")

  missing_variables <- setdiff(ordinary_variables, names(data))

  if (length(missing_variables) > 0L) {
    stop(
      "The following variables in `capture_formula` are not present in `data`: ",
      paste0("`", missing_variables, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  term_components <- unique(unlist(
    strsplit(term_labels, ":", fixed = TRUE)
  ))

  allowed_components <- formula_variables

  invalid_components <- setdiff(
    term_components,
    allowed_components
  )

  if (length(invalid_components) > 0L) {
    stop(
      "`capture_formula` may contain only untransformed variable names ",
      "and `.latent`. Invalid terms: ",
      paste0("`", invalid_components, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  main_effects <- term_labels[!grepl(":", term_labels, fixed = TRUE)]

  missing_capture_effects <- setdiff(capture_names, main_effects)

  if (length(missing_capture_effects) > 0L) {
    stop(
      "The following capture variables are missing as main effects in ",
      "`capture_formula`: ",
      paste0("`", missing_capture_effects, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  uses_latent <- ".latent" %in% formula_variables

  if (latent_classes == 1L && uses_latent) {
    stop(
      "`.latent` cannot be used when `latent_classes = 1`.",
      call. = FALSE
    )
  }

  if (
    latent_classes > 1L &&
      !".latent" %in% main_effects
  ) {
    stop(
      "`.latent` must be included as a main effect when ",
      "`latent_classes` is greater than 1.",
      call. = FALSE
    )
  }

  interaction_terms <- term_labels[
    grepl(":", term_labels, fixed = TRUE)
  ]

  for (term in interaction_terms) {
    components <- strsplit(
      term,
      split = ":",
      fixed = TRUE
    )[[1L]]

    missing_component_effects <- setdiff(
      components,
      main_effects
    )

    if (length(missing_component_effects) > 0L) {
      stop(
        "All components of an interaction must also appear as main effects. ",
        "The term `",
        term,
        "` is missing the following main effects: ",
        paste0(
          "`",
          missing_component_effects,
          "`",
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    capture_components <- intersect(
      components,
      capture_names
    )

    if (length(capture_components) == 0L) {
      next
    }

    non_capture_components <- setdiff(
      components,
      capture_names
    )

    source_only_interaction <- length(non_capture_components) == 0L

    source_latent_interaction <- length(capture_components) == 1L &&
      identical(non_capture_components, ".latent")

    if (
      !source_only_interaction &&
        !source_latent_interaction
    ) {
      stop(
        "Capture variables may interact only with other capture variables ",
        "or individually with `.latent`. Invalid term: `",
        term,
        "`.",
        call. = FALSE
      )
    }

  }

  structure(
    list(
      formula = capture_formula,
      terms = terms_object,
      term_labels = term_labels,
      variables = ordinary_variables,
      capture_names = capture_names,
      uses_latent = uses_latent
    ),
    class = "crc_capture_formula"
  )

}