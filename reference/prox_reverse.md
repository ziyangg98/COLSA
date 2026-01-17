# Proximal Operator (Reverse)

Constructs a transformation matrix to map the coefficients of Bernstein
polynomial basis functions from one dimension to another, ensuring
consistency with the properties of Bernstein polynomial bases. Note that
the original dimension must be greater than or equal to the target
dimension.

## Usage

``` r
prox_reverse(q, p)
```

## Arguments

- q:

  Integer. The number of basis functions in the original dimension.

- p:

  Integer. The number of basis functions in the target dimension.

## Value

A matrix representing the transformation from the original dimension to
the target dimension.

## Details

This function iteratively constructs the transformation matrix by
expanding the basis functions from `p` to `q`, while preserving the
mathematical properties of Bernstein polynomials.

## See also

[`prox_forward`](http://gongziyang.com/COLSA/reference/prox_forward.md)
