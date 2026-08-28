# Fit a Capture-Recapture Model

Fit a Capture-Recapture Model

## Usage

``` r
crc_fit(
  data,
  captures,
  capture_formula,
  outcome = NULL,
  outcome_formula = NULL,
  outcome_dist = c("ztnegbin", "ztpois"),
  misclass = NULL,
  latent_classes = 1L,
  control = NULL,
  verbose = FALSE
)
```

## Arguments

- data:

  A `data.frame` or `data.table` with one row per observed unit.

- captures:

  A one-sided formula identifying binary capture indicators, e.g.,
  `~ source_1 + source_2 + source_3`.

- capture_formula:

  A one-sided formula specifying the log-linear capture model. The
  special term `.latent` denotes the latent class.

- outcome:

  Either `NULL`, a single column name for an exactly observed outcome,
  or a named character vector
  `c(lower = "lower_column", upper = "upper_column")`.

- outcome_formula:

  A one-sided formula for the mean of the outcome distribution. Must be
  `NULL` when `outcome = NULL`.

- outcome_dist:

  The zero-truncated outcome distribution. Currently, `"ztnegbin"` and
  `"ztpois"` are supported.

- misclass:

  A specification of misclassification mechanisms.

- latent_classes:

  A positive integer giving the number of latent classes.

- control:

  `NULL` or a named list controlling initialization and model fitting.
  If `NULL`, the default values are used. Supported elements are:

  - `init_alpha`: A finite positive number controlling the concentration
    of the symmetric Dirichlet distribution used to initialize
    latent-class weights. The default is `20`.

  - `em`: A named list with `max_iter`, the maximum number of EM
    iterations (default `1000L`), and `tolerance`, the convergence
    tolerance (default `1e-6`).

  - `capture`: A named list with `max_iter`, the maximum number of
    iterations used in the capture-model M-step (default `100L`), and
    `tolerance`, its convergence tolerance (default `1e-8`).

  - `outcome`: A named list with `max_iter`, the maximum number of
    iterations used in the outcome-model M-step (default `1000L`), and
    `relative_tolerance`, its relative tolerance (default `1e-8`).

  Unknown control elements are rejected.

- verbose:

  A logical value indicating whether progress information should be
  displayed during model fitting. The default is `FALSE`. Currently
  unsupported.

## Value

An object of class `"crc_fit"`.
