# Proximal Operator (Forward)

Constructs a transformation matrix to map the coefficients of Bernstein
polynomial basis functions from one dimension to another, ensuring
consistency with the properties of Bernstein polynomial bases. Note that
the original dimension must be less than or equal to the target
dimension.

## Usage

``` r
prox_forward(p, q)
```

## Arguments

- p:

  Integer. The number of basis functions in the original dimension.

- q:

  Integer. The number of basis functions in the target dimension.

## Value

A matrix representing the transformation from the original dimension to
the target dimension.

## Details

This function iteratively constructs the transformation matrix by
expanding the basis functions from `p` to `q`, while preserving the
mathematical properties of Bernstein polynomials.

## See also

[`prox_reverse`](http://gongziyang.com/COLSA/reference/prox_reverse.md)
