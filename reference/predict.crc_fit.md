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
