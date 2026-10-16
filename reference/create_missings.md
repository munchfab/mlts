# Create Missings for Approximation of Continuous Time Dynamic Models

Create Missings for Approximation of Continuous Time Dynamic Models

## Usage

``` r
create_missings(data, tinterval, id, time, btw_vars = NULL)
```

## Arguments

- data:

  An object of class `data.frame` (or one that can be coerced to that
  class) containing data of all variables used in the model.

- tinterval:

  The step interval for approximation for a continuous time DSEM. The
  smaller the step interval, the better the approximation.

- id:

  The variable in `data` that identifies the person or observational
  unit (as character).

- time:

  The variable in `data` that contains the (continuous) time (as
  string).

- btw_vars:

  The names of between-level variables in the data to be added in newly
  created rows with NAs.

## Value

A `data.frame` with missings imputed for use in
[`mlts_fit`](https://munchfab.github.io/mlts/reference/mlts_fit.md).

## Examples

``` r
# create some data for example
data <- data.frame(
  id = rep(c(1, 2), each = 4),
  time = c(0, 3, 4, 6,
           1, 4, 5, 7)
)

# create missings to approximate continuous time process
create_missings(
  data = data, id = "id", time = "time",
  tinterval = 1 # use time interval of 1 minute
)
#>    num_id time int_time id
#> 1       1    0        1  1
#> 2       1   NA        2  1
#> 3       1   NA        3  1
#> 4       1    3        4  1
#> 5       1    4        5  1
#> 6       1   NA        6  1
#> 7       1    6        7  1
#> 8       2    1        1  2
#> 9       2   NA        2  2
#> 10      2   NA        3  2
#> 11      2    4        4  2
#> 12      2    5        5  2
#> 13      2   NA        6  2
#> 14      2    7        7  2
```
