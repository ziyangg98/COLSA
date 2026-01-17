# Calculate the Akaike Information Criterion (AIC) for a 'colsa' object

This function computes the AIC for a given 'colsa' object. The AIC is a
measure of the relative quality of a statistical model for a given
dataset. It is calculated as: \\-2 \* logLik + 2 \* k\\, where `k` is
the number of parameters in the model.

## Usage

``` r
# S3 method for class 'colsa'
AIC(object, ...)
```

## Arguments

- object:

  An object of class `"colsa"` for which the AIC is to be calculated.

- ...:

  Additional arguments (currently not used).

## Value

A numeric value representing the AIC of the model.

## See also

[`logLik.colsa`](http://gongziyang.com/COLSA/reference/logLik.colsa.md),
[`BIC.colsa`](http://gongziyang.com/COLSA/reference/BIC.colsa.md)
