# Simulated Survival Data with Mixed Covariates

This dataset was generated for a simulation study to evaluate survival
models incorporating both continuous and categorical covariates.

## Usage

``` r
data(sim)
```

## Format

A data frame with 6,000 rows and 10 variables:

- id:

  Integer identifier for each subject.

- time:

  Observed event or censoring time.

- status:

  Event indicator: 1 for event, 0 for right-censored.

- x1:

  First continuous covariate, generated from a bivariate normal
  distribution.

- x2:

  Second continuous covariate, also from the bivariate normal
  distribution.

- x31:

  First binary categorical covariate, generated from a Bernoulli
  distribution.

- x42:

  Dummy variable for level 2 of the second categorical covariate
  (reference is level 1).

- x43:

  Dummy variable for level 3 of the second categorical covariate.

- x44:

  Dummy variable for level 4 of the second categorical covariate.

- group:

  An integer (1 to 6) indicating the dataset partition group.

## Source

Simulation study by the package authors.

## Details

The two continuous covariates `x1` and `x2` were independently drawn
from a bivariate normal distribution. The binary covariate `x31` was
generated from a Bernoulli distribution. The second categorical
covariate had four levels and was drawn from a multinomial distribution,
conditional on `x31`. The final covariates `x42`–`x44` are dummy
variables representing levels 2–4 (level 1 is the reference).

Event times were generated from a mixture of two Weibull distributions
with shape parameters 3 and 5, and scale parameters 10 and 20,
respectively. Right censoring was imposed using censoring times drawn
from an exponential distribution with rate 3.

The true regression coefficients were: \$\$\boldsymbol{\beta} = (0.15,
-0.15, 0.3, 0.3, 0.3, 0.3)^\top\$\$

The complete dataset includes six subsets: the first three contain 1,500
observations each, and the remaining three contain 500 each. These
subsets are indicated by the `group` variable.

## Examples

``` r
data(sim)
head(sim)
#>   id       time status        x1       x2 x31 x42 x43 x44 group
#> 1  1 0.09887412      0 10.958910 5.938538   1   0   0   1     1
#> 2  2 0.05130399      0  7.499272 6.305020   1   1   0   0     1
#> 3  3 0.26820027      1  4.589824 4.250194   1   0   0   1     1
#> 4  4 0.14401458      1  4.528155 5.301456   1   0   1   0     1
#> 5  5 0.24992644      0  5.469901 5.624268   0   0   0   1     1
#> 6  6 0.40221177      1  5.063252 6.157249   1   0   0   1     1
```
