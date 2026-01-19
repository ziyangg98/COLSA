# Compute Integrated Baseline Hazard and Its Derivatives

Computes the integrated baseline hazard and its derivatives for multiple
time points using Gaussian quadrature and Bernstein polynomial basis.

## Usage

``` r
int_basehaz(times, alpha, n_basis, boundary, n_nodes = 5)
```

## Arguments

- times:

  Numeric vector of time points (upper limits of integration).

- alpha:

  Numeric vector of baseline hazard coefficients.

- n_basis:

  Integer number of Bernstein polynomial basis functions.

- boundary:

  Numeric vector of length 2 for boundary knots.

- n_nodes:

  Integer number of quadrature nodes (default 5).

## Value

A list containing:

- cbh:

  Numeric vector of cumulative baseline hazards.

- dcbh:

  Matrix of gradients (n_times x n_basis).

- d2cbh:

  Matrix of vectorized Hessians (n_times x n_basis^2).
