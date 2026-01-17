# Collaborative Inference for Cox Proportional Hazards Model in Distributed Data

Performs collaborative inference for the Cox proportional hazards model,
tailored for distributed data environments.

## Usage

``` r
colsa(formula, data, n_basis, boundary, scale = 2, init = "zero", ...)
```

## Arguments

- formula:

  A formula object specifying the model. The response must be a survival
  object created using the `Surv` function from the `survival` package.

- data:

  A data frame containing the variables in the model.

- n_basis:

  An integer specifying the number of basis functions to use.

- boundary:

  A numeric vector specifying the boundary for the basis functions. It
  should be a two-element vector with the lower and upper bounds.

- scale:

  A numeric value specifying the scaling factor for the number of
  pre-estimation basis functions. Default is 2.0.

- init:

  A character string specifying the initialization method. `"zero"`
  (default) uses zero initialization for fast computation. `"flexsurv"`
  uses `flexsurvspline` for better initial values but slower
  computation.

- ...:

  Additional arguments passed to the optimization function.

## Value

An object of class `"colsa"` containing the following components:

- logLik:

  The log-likelihood of the fitted model.

- theta:

  The estimated model parameters.

- hessian:

  The Hessian of the objective function at the solution.

- n_basis:

  The number of basis functions used.

- n_features:

  The number of features in the model.

- n_samples:

  The number of samples in the dataset.

- formula:

  The formula used to fit the model.

- boundary:

  The boundary for the basis functions.

- scale:

  The scaling factor for the number of pre-estimation basis functions.

- call:

  The matched call.

## Details

This function employs the `nlm` optimization method to estimate the
model parameters. If the optimization fails to converge, an error is
raised. During the pre-estimation stage, parameters are projected to
mitigate bias introduced in the early stage.

## See also

[`update.colsa`](http://gongziyang.com/COLSA/reference/update.colsa.md),
[`summary.colsa`](http://gongziyang.com/COLSA/reference/summary.colsa.md).

## Examples

``` r
formula <- Surv(time, status) ~ x1 + x2 + x31 + x42 + x43 + x44
boundary <- c(0, max(sim$time))
df_sub <- sim[sim$group == 1, , drop = FALSE]
fit <- colsa(formula, df_sub, n_basis = 3, boundary = boundary)
for (batch in 2:6) {
  df_sub <- sim[sim$group == batch, , drop = FALSE]
  fit <- update(fit, df_sub, n_basis = "auto")
}
summary(fit)
#> Call:
#> update.colsa(object = fit, newdata = df_sub, n_basis = "auto")
#> 
#> Number of basis functions:  6 
#> 
#>         coef exp(coef)       se      z        p    
#> x1   0.15677   1.16973  0.01016 15.430  < 2e-16 ***
#> x2  -0.17309   0.84106  0.02232 -7.756 8.74e-15 ***
#> x31  0.38400   1.46815  0.05779  6.645 3.04e-11 ***
#> x42  0.27318   1.31413  0.09090  3.005 0.002654 ** 
#> x43  0.29061   1.33725  0.08325  3.491 0.000482 ***
#> x44  0.17502   1.19128  0.08217  2.130 0.033167 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#>     exp(coef) exp(-coef) lower .95 upper .95
#> x1  1.1697    0.8549     1.1467    1.1933   
#> x2  0.8411    1.1890     0.8051    0.8787   
#> x31 1.4681    0.6811     1.3109    1.6442   
#> x42 1.3141    0.7610     1.0997    1.5704   
#> x43 1.3372    0.7478     1.1359    1.5743   
#> x44 1.1913    0.8394     1.0141    1.3994   
```
