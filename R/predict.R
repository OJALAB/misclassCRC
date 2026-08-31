#' Predict Population Size and Outcome Total
#' 
#' Computes estimates of the population size and,
#' when an outcome model was fitted, the total outcome.
#' 
#' @param object An object of class `"crc_fit"`.
#' @param by `NULL` or a one-sided formula specifying grouping variables.
#' Grouping variables must be non-capture categorical variables included in the
#' capture model. `.latent` may be used when the model has multiple
#' latent classes.
#' @param ... Additional arguments. Currently unused.
#' 
#' @return
#' An object of class `"crc_pred"`, represented as a list containing:
#' 
#' - `population_size` -- the estimated population size,
#' - `outcome_total` -- the estimated outcome total (or `NULL`
#' when no outcome model was fitted),
#' - `by_groups` -- a `data.table` containing grouped estimates
#' (or `NULL` when `by = NULL`).
#' 
#' @import data.table
#' 
#' @export
predict.crc_fit <- function(
  object,
  by = NULL,
  ...
) {

  if (!inherits(object, "crc_fit") || !isTRUE(object$fitted)) {
    stop("`object` must be a fitted `crc_fit` object.", call. = FALSE)
  }

  groups <- if (is.null(by)) {
    character()
  } else {
    if (!inherits(by, "formula") || length(by) != 2L) {
      stop("`by` must be `NULL` or a one-sided formula.", call. = FALSE)
    }

    by_terms <- terms(by)
    groups <- all.vars(by)

    if (!identical(attr(by_terms, "term.labels"), groups)) {
      stop(
        "`by` may contain only untransformed variable names joined by `+`.",
        call. = FALSE
      )
    }

    groups
  }

  if (!is.null(by) && !length(groups)) {
    stop("`by` must identify at least one grouping variable.", call. = FALSE)
  }

  allowed <- setdiff(object$capture_model$variables, object$captures$names)
  if (object$latent_classes > 1L) {
    allowed <- c(allowed, ".latent")
  }
  invalid <- setdiff(groups, allowed)
  if (length(invalid)) {
    stop(
      "Invalid grouping variables: ",
      paste(invalid, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  matrices <- object$model_matrices
  observed <- copy(matrices$states$data[, groups, with = FALSE])
  observed[, "population_size" := object$em_fit$state_weights]

  outcome_fit <- object$em_fit$outcome_fit
  if (!is.null(outcome_fit)) {
    sigma <- outcome_fit$sigma
    log_tail <- function(q, mu, shifted = FALSE) {
      if (object$outcome_dist == "ztpois") {
        return(ppois(q, lambda = mu, lower.tail = FALSE, log.p = TRUE))
      }
      size <- 1 / sigma
      if (shifted) {
        mu <- mu * (size + 1) / size
        size <- size + 1
      }
      pnbinom(q, mu = mu, size = size, lower.tail = FALSE, log.p = TRUE)
    }
    conditional_mean <- function(outcome, mu) {
      exact <- outcome$status == "exact"
      right <- outcome$status == "right_censored"
      interval <- outcome$status == "interval_censored"
      value <- exp(log(mu) - log_tail(0, mu))
      value[exact] <- outcome$lower[exact]
      value[right] <- exp(
        log(mu[right]) +
          log_tail(outcome$lower[right] - 2, mu[right], TRUE) -
          log_tail(outcome$lower[right] - 1, mu[right])
      )
      value[interval] <- exp(
        log(mu[interval]) +
          log_sub_crc(
            log_tail(outcome$lower[interval] - 2, mu[interval], TRUE),
            log_tail(outcome$upper[interval] - 1, mu[interval], TRUE)
          ) -
          log_sub_crc(
            log_tail(outcome$lower[interval] - 1, mu[interval]),
            log_tail(outcome$upper[interval], mu[interval])
          )
      )
      value
    }
    observed[,
      "outcome_total" := get("population_size") *
        conditional_mean(
          matrices$states$outcome,
          outcome_fit$fitted_mean
        )
    ]
  }

  unseen_cells <- matrices$capture$unobserved
  unseen <- copy(matrices$capture$cells[unseen_cells, groups, with = FALSE])
  unseen[,
    "population_size" := object$em_fit$capture_fit$fitted_mean[unseen_cells]
  ]
  if (!is.null(outcome_fit)) {
    x <- model.matrix(
      object$outcome_model$terms,
      matrices$capture$cells[unseen_cells]
    )
    mu <- exp(drop(x %*% outcome_fit$coefficients))
    if (anyNA(mu) || any(!is.finite(mu)) || any(mu <= 0)) {
      stop(
        "The outcome model produced invalid means for unobserved cells.",
        call. = FALSE
      )
    }
    unseen[,
      "outcome_total" := get("population_size") * exp(log(mu) - log_tail(0, mu))
    ]
  }

  result <- rbindlist(list(observed, unseen), use.names = TRUE)
  measures <- intersect(c("population_size", "outcome_total"), names(result))
  by_groups <- if (!is.null(by)) {
    result[, lapply(.SD, sum), by = groups, .SDcols = measures]
  } else NULL

  structure(
    list(
      population_size = sum(result[["population_size"]]),
      outcome_total = if (is.null(outcome_fit)) NULL else sum(result[["outcome_total"]]),
      by_groups = by_groups
    ),
    class = "crc_pred"
  )

}