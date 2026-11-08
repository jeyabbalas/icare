# Print and plot model validation results

[`print()`](https://rdrr.io/r/base/print.html) reports the study,
calibration in the large (E/O), the Hosmer-Lemeshow and relative-risk
tests, the AUC and the results of each risk-score group of a
[`ModelValidation()`](https://jeyabbalas.github.io/icare/reference/ModelValidation.md)
result. [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws
them with
[`plotModelValidation()`](https://jeyabbalas.github.io/icare/reference/plotModelValidation.md).

## Usage

``` r
# S3 method for class 'icareValid'
print(x, digits = 3, ...)

# S3 method for class 'icareValid'
plot(x, ...)
```

## Arguments

- x:

  Result of
  [`ModelValidation()`](https://jeyabbalas.github.io/icare/reference/ModelValidation.md).

- digits:

  Number of decimals of ratios and the AUC.

- ...:

  Passed to
  [`plotModelValidation()`](https://jeyabbalas.github.io/icare/reference/plotModelValidation.md)
  by [`plot()`](https://rdrr.io/r/graphics/plot.default.html); ignored
  by [`print()`](https://rdrr.io/r/base/print.html).

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns the
figure as
[`plotModelValidation()`](https://jeyabbalas.github.io/icare/reference/plotModelValidation.md)
does.
