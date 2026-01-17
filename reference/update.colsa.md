# Update a COLSA Model with New Data

This function updates a COLSA model with new data, allowing for dynamic
adjustment of the model parameters and basis functions.

## Usage

``` r
# S3 method for class 'colsa'
update(object, newdata, n_basis = "auto", alpha = 1, nu = 0.2, ...)
```

## Arguments

- object:

  An object of class `"colsa"` representing the fitted COLSA model.

- newdata:

  A `data.frame` containing the new data to update the model. It must
  include the variables specified in the original model formula.

- n_basis:

  Either a numeric value specifying the number of basis functions or
  `"auto"` to automatically determine the number of basis functions
  based on the sample size and scaling parameters. Defaults to `"auto"`,
  which means the number of basis functions will be determined by
  n_basis = round(alpha \* n_samples^nu).

- alpha:

  A positive numeric value controlling the growth rate of the number of
  basis functions when `n_basis = "auto"`. Defaults to 1.0.

- nu:

  A positive numeric value controlling the scaling exponent for
  determining the number of basis functions when `n_basis = "auto"`.
  Defaults to 0.2.

- ...:

  Additional arguments passed to the optimization function.

## Value

An updated object of class `"colsa"` with the following components:

- `logLik`: The log-likelihood of the updated model.

- `theta`: The updated parameter estimates.

- `hessian`: The updated Hessian matrix.

- `n_samples`: The updated total number of samples.

- `n_basis`: A vector tracking the number of basis functions used at
  each update step.

- `call`: The matched call of the update function.

## Details

The function updates the COLSA model by incorporating new data and
adjusting the number of basis functions dynamically. It ensures that the
model parameters and Hessian matrix are updated accordingly. The
optimization process is performed using the `nlm` method to maximize the
log-likelihood of the updated model.
