# Get Standardized Estimates for an mlts Model

Get Standardized Estimates for an mlts Model

## Usage

``` r
mlts_standardized(
  object,
  what = c("between", "within", "both"),
  digits = 3,
  prob = 0.95,
  add_cluster_std = FALSE,
  get_samples = FALSE
)
```

## Arguments

- object:

  `mltsfit`. Output of
  [`mlts_model`](https://munchfab.github.io/mlts/reference/mlts_model.md)
  and related functions.

- what:

  character. Get between-level standardized estimates
  (`what = "between"`, the default), within-level standardized estimates
  averaged over clusters (`what = "within"`), or both (`what = "both"`).

- digits:

  Number of digits. Default is 3.

- prob:

  A value between 0 and 1 to indicate the width of the credible
  interval. Default is .95.

- add_cluster_std:

  logical. If `what = "within"`, within-level standardized effects for
  each cluster are included in the output (defaults to `FALSE`).

- get_samples:

  logical. For internal use.

## Value

A `list` containing between- and within-level standardized parameters.

## Examples

``` r
# \donttest{
# build simple vector-autoregressive mlts model for two time-series variables
var_model <- mlts_model(q = 2)

# fit model with (artificial) dataset ts_data
fit <- mlts_fit(
  model = var_model,
  data = ts_data,
  ts = c("Y1", "Y2"), # time-series variables
  id = "ID", # identifier variable
  time = "time", # time variable
  tinterval = 1, # interval for approximation of continuous-time dynamic model,
  monitor_person_pars = TRUE # person parameters need to be sampled for standardization
)
#> Time series variables as indicated by parameter subscripts: 
#>    1 --> Y1
#>    2 --> Y2
#> 
#> SAMPLING FOR MODEL 'VAR_manifestCovsFix' NOW (CHAIN 1).
#> 
#> SAMPLING FOR MODEL 'VAR_manifestCovsFix' NOW (CHAIN 2).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.001387 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 13.87 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.001332 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 13.32 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 1: Iteration:   1 / 500 [  0%]  (Warmup)
#> Chain 2: Iteration:   1 / 500 [  0%]  (Warmup)
#> Chain 2: Iteration:  50 / 500 [ 10%]  (Warmup)
#> Chain 1: Iteration:  50 / 500 [ 10%]  (Warmup)
#> Chain 1: Iteration: 100 / 500 [ 20%]  (Warmup)
#> Chain 2: Iteration: 100 / 500 [ 20%]  (Warmup)
#> Chain 1: Iteration: 150 / 500 [ 30%]  (Warmup)
#> Chain 2: Iteration: 150 / 500 [ 30%]  (Warmup)
#> Chain 2: Iteration: 200 / 500 [ 40%]  (Warmup)
#> Chain 1: Iteration: 200 / 500 [ 40%]  (Warmup)
#> Chain 2: Iteration: 250 / 500 [ 50%]  (Warmup)
#> Chain 2: Iteration: 251 / 500 [ 50%]  (Sampling)
#> Chain 2: Iteration: 300 / 500 [ 60%]  (Sampling)
#> Chain 1: Iteration: 250 / 500 [ 50%]  (Warmup)
#> Chain 1: Iteration: 251 / 500 [ 50%]  (Sampling)
#> Chain 2: Iteration: 350 / 500 [ 70%]  (Sampling)
#> Chain 2: Iteration: 400 / 500 [ 80%]  (Sampling)
#> Chain 1: Iteration: 300 / 500 [ 60%]  (Sampling)
#> Chain 2: Iteration: 450 / 500 [ 90%]  (Sampling)
#> Chain 1: Iteration: 350 / 500 [ 70%]  (Sampling)
#> Chain 2: Iteration: 500 / 500 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 25.597 seconds (Warm-up)
#> Chain 2:                8.949 seconds (Sampling)
#> Chain 2:                34.546 seconds (Total)
#> Chain 2: 
#> Chain 1: Iteration: 400 / 500 [ 80%]  (Sampling)
#> Chain 1: Iteration: 450 / 500 [ 90%]  (Sampling)
#> Chain 1: Iteration: 500 / 500 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 27.792 seconds (Warm-up)
#> Chain 1:                8.122 seconds (Sampling)
#> Chain 1:                35.914 seconds (Total)
#> Chain 1: 
#> Warning: There were 151 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

# inspect standardized parameter estimates
mlts_standardized(fit)
#> $`Between-level standardized`
#> data frame with 0 columns and 0 rows
#> 
# }
```
