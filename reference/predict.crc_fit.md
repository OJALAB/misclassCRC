# Predict Population Size and Outcome Total

Computes estimates of the population size and, when an outcome model was
fitted, the total outcome.

## Usage

``` r
# S3 method for class 'crc_fit'
predict(object, by = NULL, ...)
```

## Arguments

- object:

  An object of class `"crc_fit"`.

- by:

  `NULL` or a one-sided formula specifying grouping variables. Grouping
  variables must be non-capture categorical variables included in the
  capture model. `.latent` may be used when the model has multiple
  latent classes.

- ...:

  Additional arguments. Currently unused.

## Value

An object of class `"crc_pred"`, represented as a list containing:

- `population_size` – the estimated population size,

- `outcome_total` – the estimated outcome total (or `NULL` when no
  outcome model was fitted),

- `by_groups` – a `data.table` containing grouped estimates (or `NULL`
  when `by = NULL`).

## See also

Other fitting functions:
[`crc_fit()`](https://ojalab.ue.poznan.pl/misclassCRC/reference/crc_fit.md),
[`print.crc_fit()`](https://ojalab.ue.poznan.pl/misclassCRC/reference/print.crc_fit.md),
[`print.summary.crc_fit()`](https://ojalab.ue.poznan.pl/misclassCRC/reference/print.summary.crc_fit.md),
[`summary.crc_fit()`](https://ojalab.ue.poznan.pl/misclassCRC/reference/summary.crc_fit.md)

Other prediction functions:
[`print.crc_pred()`](https://ojalab.ue.poznan.pl/misclassCRC/reference/print.crc_pred.md),
[`print.summary.crc_pred()`](https://ojalab.ue.poznan.pl/misclassCRC/reference/print.summary.crc_pred.md),
[`summary.crc_pred()`](https://ojalab.ue.poznan.pl/misclassCRC/reference/summary.crc_pred.md)
