# Print Capture-Recapture Predictions

Displays the population-size and outcome-total estimates computed by
[`predict.crc_fit()`](https://ojalab.ue.poznan.pl/misclassCRC/reference/predict.crc_fit.md),
including grouped estimates when requested.

## Usage

``` r
# S3 method for class 'crc_pred'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"crc_pred"`.

- digits:

  The number of significant digits to display.

- ...:

  Additional arguments. Currently unused.

## Value

`x`, invisibly.

## See also

Other prediction functions:
[`predict.crc_fit()`](https://ojalab.ue.poznan.pl/misclassCRC/reference/predict.crc_fit.md),
[`print.summary.crc_pred()`](https://ojalab.ue.poznan.pl/misclassCRC/reference/print.summary.crc_pred.md),
[`summary.crc_pred()`](https://ojalab.ue.poznan.pl/misclassCRC/reference/summary.crc_pred.md)
