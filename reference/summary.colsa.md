# Summarize a 'colsa' Object

This function provides a summary for objects of class `"colsa"`.

## Usage

``` r
# S3 method for class 'colsa'
summary(object, conf.int = 0.95, ...)
```

## Arguments

- object:

  An object of class `"colsa"`.

- conf.int:

  A numeric value between 0 and 1 specifying the confidence level for
  the confidence intervals. Default is 0.95.

- ...:

  Additional arguments (currently not used).

## Value

A list of class `"summary.colsa"` containing the following components:

- call:

  The function call that created the `"colsa"` object.

- coefficients:

  A matrix with columns for the coefficients, exponentiated
  coefficients, standard errors, z-scores, and p-values.

- conf.int:

  A matrix with columns for the exponentiated coefficients, their
  reciprocals, and the lower and upper confidence interval bounds.

- n_basis:

  The number of basis functions used in the model.

- boundary:

  The boundary information from the `"colsa"` object.

## Details

The function calculates the z-scores and p-values for the coefficients,
as well as confidence intervals for the exponentiated coefficients.
