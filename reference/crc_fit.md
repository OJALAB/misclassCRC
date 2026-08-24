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
  control = NULL
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

  An optional control list; only `init_alpha` is currently supported.

## Value

An object of class `"crcfit"`.
