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

  categorical_variables <- setdiff(
    ordinary_variables,
    capture_names
  )

  validate_categorical_variables(
    variables = categorical_variables,
    data = data,
    context = "capture_formula"
  )

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

#' @importFrom stats terms
#' 
#' @noRd
parse_outcome_formula <- function(
  outcome_formula,
  data,
  outcome,
  latent_classes
) {

  if (is.null(outcome)) {
    if (!is.null(outcome_formula)) {
      stop(
        "`outcome_formula` must be `NULL` when `outcome` is `NULL`.",
        call. = FALSE
      )
    }
    return(NULL)
  }

  if (
    !inherits(outcome_formula, "formula") ||

      length(outcome_formula) != 2L
  ) {
    stop(
      "`outcome_formula` must be a one-sided formula.",
      call. = FALSE
    )
  }

  terms_object <- terms(outcome_formula)
  term_labels <- attr(terms_object, "term.labels")
  
  if (attr(terms_object, "intercept") != 1L) {
    stop(
      "`outcome_formula` must include an intercept.",
      call. = FALSE
    )
  }
  
  formula_variables <- all.vars(outcome_formula)
  ordinary_variables <- setdiff(formula_variables, ".latent")

  missing_variables <- setdiff(
    ordinary_variables,
    names(data)
  )

  if (length(missing_variables) > 0L) {
    stop(
      "The following variables in `outcome_formula` are not present in `data`: ",
      paste0("`", missing_variables, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  validate_categorical_variables(
    variables = ordinary_variables,
    data = data,
    context = "outcome_formula"
  )

  term_components <- unique(unlist(
    strsplit(term_labels, ":", fixed = TRUE)
  ))

  invalid_components <- setdiff(
    term_components,
    formula_variables
  )

  if (length(invalid_components) > 0L) {
    stop(
      "`outcome_formula` may contain only untransformed variable names ",
      "and `.latent`. Invalid terms: ",
      paste0("`", invalid_components, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  uses_latent <- ".latent" %in% formula_variables

  if (latent_classes == 1L && uses_latent) {
    stop(
      "`.latent` cannot be used in `outcome_formula` when ",
      "`latent_classes = 1`.",
      call. = FALSE
    )
  }

  main_effects <- term_labels[
    !grepl(":", term_labels, fixed = TRUE)
  ]

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
  }

  structure(
    list(
      formula = outcome_formula,
      terms = terms_object,
      term_labels = term_labels,
      variables = ordinary_variables,
      uses_latent = uses_latent
    ),
    class = "crc_outcome_formula"
  )
  
}

#' @import data.table
#' @importFrom stats terms
#' 
#' @noRd
parse_misclass <- function(
  misclass,
  data,
  capture_names
) {

  if (is.null(misclass)) {
    return(NULL)
  }

  if (
    !is.list(misclass) ||
      is.null(names(misclass)) ||
      anyNA(names(misclass)) ||
      any(!nzchar(names(misclass))) ||
      anyDuplicated(names(misclass)) > 0L
  ) {
    stop(
      "`misclass` must be a named list of misclassification specifications.",
      call. = FALSE
    )
  }

  variables <- names(misclass)

  missing_variables <- setdiff(
    variables,
    names(data)
  )

  if (length(missing_variables) > 0L) {
    stop(
      "The following variables specified in `misclass` are not present in `data`: ",
      paste0("`", missing_variables, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  parsed <- vector(
    mode = "list",
    length = length(misclass)
  )

  names(parsed) <- variables

  for (variable in variables) {

    spec <- misclass[[variable]]

    if (!is.list(spec)) {
      stop(
        "The misclassification specification for `",
        variable,
        "` must be a list.",
        call. = FALSE
      )
    }

    if (
      is.null(names(spec)) ||
        anyNA(names(spec)) ||
        any(!nzchar(names(spec))) ||
        anyDuplicated(names(spec)) > 0L
    ) {
      stop(
        "The misclassification specification for `",
        variable,
        "` must have unique, non-empty element names.",
        call. = FALSE
      )
    }

    allowed_names <- c(
      "matrix",
      "true_if"
    )

    unknown_names <- setdiff(
      names(spec),
      allowed_names
    )

    if (length(unknown_names) > 0L) {
      stop(
        "Unknown elements in the misclassification specification for `",
        variable,
        "`: ",
        paste0("`", unknown_names, "`", collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (!"matrix" %in% names(spec)) {
      stop(
        "The misclassification specification for `",
        variable,
        "` must contain `matrix`.",
        call. = FALSE
      )
    }

    Pi <- spec$matrix

    if (
      !is.matrix(Pi) ||
        !is.numeric(Pi)
    ) {
      stop(
        "The misclassification matrix for `",
        variable,
        "` must be a numeric matrix.",
        call. = FALSE
      )
    }

    if (nrow(Pi) != ncol(Pi)) {
      stop(
        "The misclassification matrix for `",
        variable,
        "` must be square.",
        call. = FALSE
      )
    }

    if (
      anyNA(Pi) ||
        any(!is.finite(Pi)) ||
        any(Pi < 0)
    ) {
      stop(
        "The misclassification matrix for `",
        variable,
        "` must contain finite non-negative values.",
        call. = FALSE
      )
    }

    if (!all(abs(rowSums(Pi) - 1) < 1e-8)) {
      stop(
        "Rows of the misclassification matrix for `",
        variable,
        "` must sum to 1.",
        call. = FALSE
      )
    }

    if (
      is.null(rownames(Pi)) ||
        is.null(colnames(Pi))
    ) {
      stop(
        "The misclassification matrix for `",
        variable,
        "` must have row and column names.",
        call. = FALSE
      )
    }

    if (
      anyNA(rownames(Pi)) ||
        anyNA(colnames(Pi)) ||
        any(!nzchar(rownames(Pi))) ||
        any(!nzchar(colnames(Pi)))
    ) {
      stop(
        "Row and column names of the misclassification matrix for `",
        variable,
        "` must be non-missing and non-empty.",
        call. = FALSE
      )
    }

    if (
      anyDuplicated(rownames(Pi)) > 0L ||
        anyDuplicated(colnames(Pi)) > 0L
    ) {
      stop(
        "Row and column names of the misclassification matrix for `",
        variable,
        "` must be unique.",
        call. = FALSE
      )
    }

    if (!identical(rownames(Pi), colnames(Pi))) {
      stop(
        "Row and column names of the misclassification matrix for `",
        variable,
        "` must contain the same category levels in the same order.",
        call. = FALSE
      )
    }

    observed_values <- data[[variable]]

    if (
      !is.factor(observed_values) &&
        !is.character(observed_values)
    ) {
      stop(
        "Misclassified variable `",
        variable,
        "` must be a factor or character variable.",
        call. = FALSE
      )
    }

    if (anyNA(observed_values)) {
      stop(
        "Misclassified variable `",
        variable,
        "` must not contain missing values.",
        call. = FALSE
      )
    }

    observed_levels <- unique(
      as.character(observed_values)
    )

    unknown_levels <- setdiff(
      observed_levels,
      colnames(Pi)
    )

    if (length(unknown_levels) > 0L) {
      stop(
        "The following observed levels of `",
        variable,
        "` are not represented in the misclassification matrix: ",
        paste0("`", unknown_levels, "`", collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    true_if <- spec$true_if

    if (is.null(true_if)) {

      true_sources <- character(0L)

      is_true <- rep(
        FALSE,
        nrow(data)
      )

    } else {

      if (
        !inherits(true_if, "formula") ||
          length(true_if) != 2L
      ) {
        stop(
          "`true_if` for `",
          variable,
          "` must be a one-sided formula.",
          call. = FALSE
        )
      }

      true_terms_object <- terms(true_if)
      true_terms <- attr(
        true_terms_object,
        "term.labels"
      )

      true_sources <- all.vars(true_if)

      if (
        !identical(true_terms, true_sources) ||

          attr(true_terms_object, "intercept") != 1L
      ) {
        stop(
          "`true_if` for `",
          variable,
          "` must contain only capture-variable names joined with `+`.",
          call. = FALSE
        )
      }

      invalid_sources <- setdiff(
        true_sources,
        capture_names
      )

      if (length(invalid_sources) > 0L) {
        stop(
          "The following variables in `true_if` for `",
          variable,
          "` are not capture variables: ",
          paste0("`", invalid_sources, "`", collapse = ", "),
          ".",
          call. = FALSE
        )
      }

      if (length(true_sources) == 0L) {
        stop(
          "`true_if` for `",
          variable,
          "` must identify at least one capture source.",
          call. = FALSE
        )
      }

      true_data <- as.data.table(data)[
        ,
        true_sources,
        with = FALSE
      ]

      is_true <- rowSums(
        as.matrix(true_data)
      ) > 0L

    }

    parsed[[variable]] <- list(
      variable = variable,
      matrix = Pi,
      levels = rownames(Pi),
      true_sources = true_sources,
      is_true = is_true
    )

  }

  structure(
    list(
      variables = variables,
      specifications = parsed
    ),
    class = "crc_misclass"
  )

}

#' @noRd
parse_control <- function(control) {

  defaults <- list(
    init_alpha = 20,
    em = list(max_iter = 1000L, tolerance = 1e-6),
    capture = list(max_iter = 100L, tolerance = 1e-8),
    outcome = list(max_iter = 1000L, relative_tolerance = 1e-8)
  )

  if (is.null(control)) {
    control <- list()
  }

  check_names <- function(x, allowed, context) {
    if (!is.list(x)) {
      stop("`", context, "` must be a list.", call. = FALSE)
    }
    if (length(x) == 0L) {
      return(invisible(TRUE))
    }
    element_names <- names(x)
    if (
      is.null(element_names) ||
        anyNA(element_names) ||
        any(!nzchar(element_names)) ||
        anyDuplicated(element_names)
    ) {
      stop("`", context, "` must have unique, non-empty names.", call. = FALSE)
    }
    unknown <- setdiff(element_names, allowed)
    if (length(unknown) > 0L) {
      stop(
        "Unknown setting in `",
        context,
        "`: `",
        paste(unknown, collapse = "`, `"),
        "`.",
        call. = FALSE
      )
    }
    invisible(TRUE)
  }

  check_names(control, names(defaults), "control")
  resolved <- defaults

  if ("init_alpha" %in% names(control)) {
    resolved$init_alpha <- control$init_alpha
  }

  for (section in c("em", "capture", "outcome")) {
    if (section %in% names(control)) {
      check_names(
        control[[section]],
        names(defaults[[section]]),
        paste0("control$", section)
      )
      resolved[[section]][names(control[[section]])] <- control[[section]]
    }
  }

  positive <- function(x, integer = FALSE) {
    is.numeric(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      is.finite(x) &&
      x > 0 &&
      (!integer || x == floor(x))
  }

  if (!positive(resolved$init_alpha)) {
    stop("`control$init_alpha` must be a finite positive number.", call. = FALSE)
  }
  if (!positive(resolved$em$max_iter, TRUE)) {
    stop("`control$em$max_iter` must be a positive integer.", call. = FALSE)
  }
  if (!positive(resolved$em$tolerance)) {
    stop("`control$em$tolerance` must be a finite positive number.", call. = FALSE)
  }
  if (!positive(resolved$capture$max_iter, TRUE)) {
    stop("`control$capture$max_iter` must be a positive integer.", call. = FALSE)
  }
  if (!positive(resolved$capture$tolerance)) {
    stop(
      "`control$capture$tolerance` must be a finite positive number.",
      call. = FALSE
    )
  }
  if (!positive(resolved$outcome$max_iter, TRUE)) {
    stop("`control$outcome$max_iter` must be a positive integer.", call. = FALSE)
  }
  if (!positive(resolved$outcome$relative_tolerance)) {
    stop(
      "`control$outcome$relative_tolerance` must be a finite positive number.",
      call. = FALSE
    )
  }

  structure(resolved, class = "crc_control")

}