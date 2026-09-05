#' Print a Capture-Recapture Model Fit
#' 
#' Displays a concise description of a model fitted by [crc_fit()],
#' including its main components and convergence status.
#' 
#' @param x An object of class `"crc_fit"`.
#' @param ... Additional arguments. Currently unused.
#' 
#' @return 
#' `x`, invisibly.
#' 
#' @export
print.crc_fit <- function(
  x,
  ...
) {

  format_integer <- function(value) {
    format(value, big.mark = ",", scientific = FALSE, trim = TRUE)
  }

  misclassified <- if (is.null(x$misclass)) {
    "none"
  } else {
    paste(x$misclass$variables, collapse = ", ")
  }

  outcome_distribution <- if (is.null(x$outcome_model)) {
    "none"
  } else {
    switch(
      x$outcome_dist,
      ztnegbin = "zero-truncated negative binomial",
      ztpois = "zero-truncated Poisson"
    )
  }

  iterations <- x$em_fit$iterations
  iteration_label <- if (iterations == 1L) "iteration" else "iterations"
  convergence <- paste(
    if (isTRUE(x$em_fit$converged)) {
      "converged after"
    } else {
      "not converged after"
    },
    format_integer(iterations),
    iteration_label
  )

  values <- c(
    "Observed units:" = format_integer(nrow(x$data)),
    "Capture sources:" = paste(x$captures$names, collapse = ", "),
    "Latent classes:" = format_integer(x$latent_classes),
    "Misclassified variables:" = misclassified,
    "Outcome distribution:" = outcome_distribution,
    "EM algorithm:" = convergence
  )

  cat("Capture-recapture model fit\n\n")
  cat(sprintf("%-25s%s\n", names(values), unname(values)), sep = "")

  invisible(x)

}

#' Print a Comprehensive Summary of a Capture-Recapture
#' Model Fit
#' 
#' @param x An object of class `"summary.crc_fit"`.
#' @param digits The number of significant digits to display.
#' @param ... Additional arguments. Currently unused.
#' 
#' @return
#' `x`, invisibly.
#' 
#' @export
print.summary.crc_fit <- function(
  x,
  digits = max(3L, getOption("digits") - 3L),
  ...
) {

  line <- function(label, value) {
    value <- format(value, digits = digits, big.mark = ",", trim = TRUE)
    cat(sprintf("  %-25s%s\n", paste0(label, ":"), value), sep = "")
  }

  cat("Capture-recapture model summary\n\nCall:\n")
  print(x$call)

  cat("\nModel structure:\n")
  for (name in names(x$model)) {
    line(name, x$model[[name]])
  }

  if (length(x$misclassification)) {
    cat("\nMisclassification:\n")
    for (name in names(x$misclassification)) {
      line(name, x$misclassification[[name]])
    }
  }

  cat("\nCapture model:\n  ")
  print(x$capture$formula)
  print(cbind(Estimate = x$capture$coefficients), digits = digits)

  if (!is.null(x$outcome)) {
    cat("\nOutcome model:\n  ")
    print(x$outcome$formula)
    print(cbind(Estimate = x$outcome$coefficients), digits = digits)
    if (!is.null(x$outcome$sigma)) line("Dispersion (sigma)", x$outcome$sigma)
  }

  cat("\nEM algorithm:\n")
  line("Status", if (isTRUE(x$converged)) "converged" else "not converged")
  line("Iterations", x$iterations)
  line("Maximum change", x$convergence$maximum)
  line("State-weight change", x$convergence$changes[["state_weights"]])
  line("Cell-count change", x$convergence$changes[["cell_counts"]])

  invisible(x)

}
