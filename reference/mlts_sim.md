# Simulate data from mlts model

Simulate data from mlts model

## Usage

``` r
mlts_sim(
  model,
  default = FALSE,
  N = NULL,
  N_G = NULL,
  TP,
  burn.in = 50,
  seed = NULL,
  seed.true = 1,
  btw.var.sds = NULL,
  exogenous = NULL
)
```

## Arguments

- model:

  `data.frame`. Output of
  [`mlts_model`](https://munchfab.github.io/mlts/reference/mlts_model.md).

- default:

  logical. If set to `TRUE`, default prior specifications are added.

- N:

  integer Number of observational units.

- N_G:

  vector of integers. Number of observational units per group.

- TP:

  integer. Number of measurements per observational unit.

- burn.in:

  integer. Length of ‘burn-in’ period.

- seed:

  integer. Seed used for data generation.

- seed.true:

  integer. Separate seed used for sampling of true population parameters
  values from plausible ranges for stationary time series.

- btw.var.sds:

  named numeric vector. Provide standard deviation(s) for all exogenous
  between-level variable(s) specified in `model`, e.g.
  (`btw.var.sds = c("covariate1" = 1)`, to set the SD of the variable
  "covariate1" to 1). Mean values of the respective variable(s) will be
  set to 0 per default.

- exogenous:

  Matrix of numeric values of exogenous variables with
  `N`\*(`TP`+`burn.in`) rows and separate columns for each variable.

## Value

An object of class `"mlts_simdata"`. The object is a list containing the
following components:

- model:

  the model object passed to `mlts_sim` with true parameter values used
  in the data generation added in the column `true.val`

- data:

  a long format `data.frame` of the generated time series data

- RE.pars:

  a `matrix` of cluster-specific true values used in the data generation

## Details

A function to generate data from an output of
[`mlts_model`](https://munchfab.github.io/mlts/reference/mlts_model.md).

## Examples

``` r
# \donttest{
# build a simple vector-autoregressive mlts model with two time-series variables
var_model <- mlts_model(q = 2)

# simulate data from this model with default true values
# (true values are randomly drawn from normal distribution)
var_data <- mlts_sim(
  model = var_model,
  N = 50, TP = 30, # number of units and number of measurements per unit
  default = TRUE # use default parameter values
)

# the data set is stored in .$data
head(var_data$data)
#>   ID time         Y1          Y2
#> 1  1    1  1.9031184 -1.97954112
#> 2  1    2  0.4627377 -2.25785898
#> 3  1    3  1.8569070 -2.23935333
#> 4  1    4  0.9839013 -0.01119517
#> 5  1    5 -0.7042657 -1.74205278
#> 6  1    6  0.3960619 -2.39583893

# individual parameter values are stored in .$RE.pars
head(var_data$RE.pars)
#>             mu_1       mu_2   phi(1)_11   phi(1)_12   phi(1)_21    phi(1)_22
#> [1,]  0.71782500 -1.3054433  0.14629806 -0.37034855 -0.23454961  0.199117316
#> [2,]  0.25497664  1.9156897 -0.07568515 -0.28399710  0.17465057  0.160338484
#> [3,]  0.44272628  2.3889042  0.15165108 -0.33393953  0.24189355  0.001282907
#> [4,] -0.03415591  0.8297450  0.33296081  0.02407496  0.35121338  0.430938860
#> [5,]  0.22458128  1.1653852 -0.07875350 -0.13737890 -0.08229474 -0.045020808
#> [6,]  0.27738081 -0.2664382  0.33690756 -0.35183497  0.09533000  0.102151359
#>      ln.sigma2_1  ln.sigma2_2
#> [1,] -0.44668021 -0.341635125
#> [2,] -0.38051325 -0.291195608
#> [3,] -0.09932916 -0.004460711
#> [4,] -0.18390499 -0.498366716
#> [5,] -0.20451911 -0.448134693
#> [6,]  0.40130511 -0.319759467

# if the mltssim-object is used in mlts_fit(), true values
# are added to the fitted object
fit <- mlts_fit(
  model = var_model,
  data = var_data,
  id = "ID", ts = c("Y1", "Y2"), time = "time"
)
#> Time series variables as indicated by parameter subscripts: 
#>    1 --> Y1
#>    2 --> Y2
#> Simulated data provided:
#> True scores used in the data generation will be added to the returned object.
#> 
#>  Note that obtaining (standardized) estimates by cluster,
#>     requires settting `monitor_person_pars = TRUE`. However, keeping the default option
#>     (`monitor_person_pars = FALSE`) can improve sampling times.
#> Warning: There were 3 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> https://mc-stan.org/misc/warnings.html#bfmi-low
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

# inspect model with true values
head(fit$pop.pars.summary)
#>    Param_stan     Param      Model  Level         Type Param_Label isRandom
#> 1 gammas[1,1]      mu_1 Structural Within Fixed effect       Trait        1
#> 2 gammas[1,2]      mu_2 Structural Within Fixed effect       Trait        1
#> 3 gammas[1,3] phi(1)_11 Structural Within Fixed effect     Dynamic        1
#> 4 gammas[1,4] phi(1)_12 Structural Within Fixed effect     Dynamic        1
#> 5 gammas[1,5] phi(1)_21 Structural Within Fixed effect     Dynamic        1
#> 6 gammas[1,6] phi(1)_22 Structural Within Fixed effect     Dynamic        1
#>   prior_type prior_location prior_scale group Par_no true.val   mean    sd
#> 1     normal              0          10     1      1     0.80  0.709 0.107
#> 2     normal              0          10     1      2     0.30  0.479 0.219
#> 3     normal              0           2     1      3     0.15  0.185 0.037
#> 4     normal              0           2     1      4    -0.20 -0.241 0.031
#> 5     normal              0           2     1      5     0.05  0.064 0.033
#> 6     normal              0           2     1      6     0.15  0.142 0.045
#>     2.5%    50%  97.5% n_eff  Rhat Bulk_ESS Tail_ESS
#> 1  0.503  0.712  0.910   582 1.012      599      324
#> 2  0.045  0.491  0.874   945 1.010     1050      289
#> 3  0.108  0.185  0.262   302 1.003      315      351
#> 4 -0.308 -0.242 -0.186   247 1.001      243      286
#> 5  0.001  0.064  0.126   349 1.005      358      342
#> 6  0.054  0.142  0.227   450 1.001      451      419
#> %]  (Sampling)
#> Chain 2: Iteration: 400 / 500 [ 80%]  (Sampling)
#> Chain 2: Iteration: 450 / 500 [ 90%]  (Sampling)
#> Chain 2: Iteration: 500 / 500 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 27.033 seconds (Warm-up)
#> Chain 2:                75.741 seconds (Sampling)
#> Chain 2:                102.774 seconds (Total)
#> Chain 2: 
# }
```
