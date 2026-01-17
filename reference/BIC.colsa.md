# Calculate Bayesian Information Criterion (BIC) for a 'colsa' Object

This function computes the Bayesian Information Criterion (BIC) for a
fitted 'colsa' model object. The BIC is a model selection criterion that
balances model fit and complexity.

## Usage

``` r
# S3 method for class 'colsa'
BIC(object, ...)
```

## Arguments

- object:

  An object of class 'colsa'. This is the fitted model for which the BIC
  is to be calculated.

- ...:

  Additional arguments (currently not used).

## Value

A numeric value representing the BIC of the model.

## Details

The BIC is calculated as: \$\$-2 \times \log(\text{Likelihood}) + k
\times \log(n)\$\$ where \\k\\ is the number of parameters in the model,
and \\n\\ is the number of samples in the dataset.

## See also

[`logLik.colsa`](http://gongziyang.com/COLSA/reference/logLik.colsa.md),
[`AIC.colsa`](http://gongziyang.com/COLSA/reference/AIC.colsa.md)
