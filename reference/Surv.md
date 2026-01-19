# Alias for Survival Function

This function creates an alias for the `Surv` function from the
`survival` package, allowing it to be used without explicitly
referencing the package namespace.

## Usage

``` r
Surv(
  time,
  time2,
  event,
  type = c("right", "left", "interval", "counting", "interval2"),
  origin = 0
)
```

## Arguments

- time:

  The follow-up time for the subject.

- time2:

  The ending time for the interval (if applicable).

- event:

  The status indicator, normally 0=alive, 1=dead. Other values are
  allowed.

- type:

  Character string specifying the type of censoring. Possible values are
  "right", "left", "interval", or "counting".

- origin:

  The origin for counting time, used when `type = "counting"`.

## Value

A `Surv` object representing the survival data.

## See also

[`Surv`](https://rdrr.io/pkg/survival/man/Surv.html)
