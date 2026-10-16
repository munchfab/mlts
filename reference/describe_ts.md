# Describe Time Series Variables by Group

Computes descriptive statistics for one or more time-series variables,
optionally within groups. For each variable (and group, if provided),
the function returns sample size, number of complete observations,
percentage of missing values, number of complete lag-1 pairs, central
tendency, dispersion, skewness, kurtosis, and floor/ceiling effects.

## Usage

``` r
describe_ts(
  data,
  ts,
  group = NULL,
  digits = 3,
  scale_min = NULL,
  scale_max = NULL,
  skew_type = 2,
  kurt_type = 2
)
```

## Arguments

- data:

  A `data.frame` containing the time-series variables.

- ts:

  A character vector of variable names (columns in `data`) representing
  the time-series to be summarized.

- group:

  Optional. A character string giving the name of a grouping variable in
  `data`. If `NULL` (default), all observations are treated as belonging
  to one group.

- digits:

  Integer indicating the number of decimal places to which numeric
  statistics should be rounded. Default is 3.

- scale_min, scale_max:

  Optional numeric values specifying the minimum and maximum of the
  measurement scale for computing floor and ceiling percentages. If
  `NULL` (default), these are taken as the minimum and maximum observed
  values in `data[[var]]` for each time-series variable.

- skew_type:

  Integer (1, 2, or 3) passed to
  [`e1071::skewness`](https://rdrr.io/pkg/e1071/man/skewness.html) that
  controls the estimator of skewness (see Joanes & Gill, 1998, and the
  `e1071` documentation). The default is `2`, which corresponds to the
  estimator used by SPSS/SAS.

- kurt_type:

  Integer (1, 2, or 3) passed to
  [`e1071::kurtosis`](https://rdrr.io/pkg/e1071/man/kurtosis.html) that
  controls the estimator of kurtosis. The default is `2`, matching the
  skewness default.

## Value

A `data.frame` with one row per variable per group, containing:

- `group` — group identifier (omitted if `group = NULL`)

- `ts_var` — name of time series variable

- `N` — total number of observations (including NAs)

- `N_comp` — number of non-missing observations

- `miss_pc` — percent missing

- `N_lag1_comp` — number of complete lag-1 pairs

- `M` — mean

- `Mdn` — median

- `SD` — standard deviation

- `Skew` — sample skewness (as returned by
  [`e1071::skewness`](https://rdrr.io/pkg/e1071/man/skewness.html) with
  the specified `skew_type`)

- `Kurtosis` — sample kurtosis (as returned by
  [`e1071::kurtosis`](https://rdrr.io/pkg/e1071/man/kurtosis.html) with
  the specified `kurt_type`)

- `floor_pc` — percent at `scale_min`

- `ceiling_pc` — percent at `scale_max`

## Details

For each time-series variable, the function removes `NA` values before
computing summary statistics. Skewness and kurtosis are obtained from
[`e1071::skewness`](https://rdrr.io/pkg/e1071/man/skewness.html) and
[`e1071::kurtosis`](https://rdrr.io/pkg/e1071/man/kurtosis.html),
respectively, using the `skew_type` and `kurt_type` arguments. When the
standard deviation is zero, missing, or when fewer than four complete
observations are available, both skewness and kurtosis are set to `NA`.

Floor and ceiling effects are defined as the percentage of non-missing
observations equal to `scale_min` and `scale_max`, respectively.

Lag-1 completeness is computed as the number of adjacent observation
pairs where both the current and next observation are non-missing. If
`group = NULL`, this computation does not account for clustering or
overnight gaps—that is, the function assumes the data form a single
continuous time series, and lag-1 pairs may span across natural breaks
in the data.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example with grouping
describe_ts(data = df, ts = c("x", "y"), group = "id")

# Example without grouping
describe_ts(data = df, ts = "x")
} # }
```
