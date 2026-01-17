# Objective Function for Cox Model Optimization

Computes the objective value, gradient, and Hessian for optimizing a Cox
proportional hazards model with spline-based baseline hazard estimation.

## Usage

``` r
objective(par, time, status, x, boundary, theta, hessian_prev)
```

## Arguments

- par:

  Numeric vector of parameters to be optimized, including spline
  coefficients for the baseline hazard (`alpha`) and regression
  coefficients (`beta`).

- time:

  Numeric vector of observed survival times.

- status:

  Numeric vector of event indicators.

- x:

  Numeric matrix of covariates (features) for the model.

- boundary:

  Numeric vector specifying the boundary knots for the splines.

- theta:

  Numeric vector of prior parameters for regularization.

- hessian_prev:

  Numeric matrix of prior Hessian for regularization.

## Value

A list containing:

- value:

  Scalar value of the objective function.

- gradient:

  Numeric vector of the gradient of the objective function.

- hessian:

  Numeric matrix of the Hessian of the objective function.

## Details

The objective function combines two components:

- A regularization term based on prior parameters (`theta`) and prior
  Hessian (`hessian`).

- A likelihood-based term for the Cox proportional hazards model, where
  the baseline hazard is modeled using Bernstein polynomials.

The function calculates the loss, gradient, and Hessian for both
components and combines them to return the overall objective value,
gradient, and Hessian. The baseline hazard is integrated using helper
functions (`int_basehaz`) to compute cumulative hazard and its
derivatives.
