# Extract Coefficients from a 'colsa' Object

This function retrieves the coefficients from a fitted 'colsa' object.

## Usage

``` r
# S3 method for class 'colsa'
coef(object, ...)
```

## Arguments

- object:

  An object of class 'colsa'. The fitted model object from which
  coefficients are to be extracted.

- ...:

  Additional arguments (currently not used).

## Value

A numeric vector containing the coefficients extracted from the 'colsa'
object.

## Details

The function verifies that the input object is of class 'colsa'. It then
calculates the coefficients by applying a transformation matrix (`prox`)
to the `theta` parameter of the object. The transformation matrix
accounts for the number of basis functions and features in the model.
