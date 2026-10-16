# Create a summary of a fitted model with class `mltsfit`

Create a summary of a fitted model with class `mltsfit`

## Usage

``` r
# S3 method for class 'mltsfit'
summary(
  object,
  priors = FALSE,
  se = FALSE,
  prob = 0.95,
  bpe = c("mean"),
  digits = 3,
  flag_signif = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class `mltsfit`.

- priors:

  Add prior information (default = FALSE).

- se:

  Logical. Should the Monte Carlo Standard Error be included in the
  summary? Defaults to `FALSE`.

- prob:

  A value between 0 and 1 to indicate the width of the credible
  interval. Default is .95.

- bpe:

  Bayesian posterior estimate can be either "mean" (the default) or the
  "median" of the posterior distribution.

- digits:

  Number of digits.

- flag_signif:

  Add significance flags based on `prob` (default = FALSE).

- ...:

  Additional arguments affecting the summary produced.

## Value

A summary of model parameters.

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
  time = "time",
  tinterval = 1 # interval for approximation of continuous-time dynamic model,
)
#> Time series variables as indicated by parameter subscripts: 
#>    1 --> Y1
#>    2 --> Y2
#> 
#>  Note that obtaining (standardized) estimates by cluster,
#>     requires settting `monitor_person_pars = TRUE`. However, keeping the default option
#>     (`monitor_person_pars = FALSE`) can improve sampling times.
#> Warning: There were 7 divergent transitions after warmup. See
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

# inspect model summary
summary(fit)
#> Time series variables as indicated by parameter subscripts: 
#>    1 --> Y1
#>    2 --> Y2
#> Data: 1100 observations in 55 IDs
#> Model convergence criteria: 
#>   Maximum Potential Scale Reduction Factor (PSR; Rhat): 1.49 (should be < 1.01)
#>   Minimum Bulk ESS: 4 (should be > 200, 100 per chain) 
#>   Minimum Tail ESS: 23 (should be > 200, 100 per chain) 
#>   Number of divergent transitions: 7 (should be 0) 
#> 
#> Posterior Summary Statistics
#> Fixed Effects:
#>                Mean    SD   2.5%  97.5%  Rhat Bulk_ESS Tail_ESS
#>         mu_1  0.777 0.154  0.449  1.075 1.003      350      199
#>         mu_2  0.341 0.180 -0.029  0.675 1.004      408      260
#>    phi(1)_11  0.355 0.039  0.281  0.435 1.035       89       84
#>    phi(1)_12 -0.107 0.045 -0.191 -0.028 1.004      298      287
#>    phi(1)_21  0.059 0.039 -0.024  0.128 1.003      123      381
#>    phi(1)_22  0.134 0.041  0.060  0.209 1.020      182      172
#>  ln.sigma2_1 -0.207 0.058 -0.330 -0.098 1.044       64      121
#>  ln.sigma2_2 -0.294 0.057 -0.405 -0.183 1.033      118      265
#>    r.zeta_12 -0.145 0.032 -0.206 -0.080 1.011      376      313
#> 
#> Random Effects SDs:
#>                     Mean    SD  2.5% 97.5%  Rhat Bulk_ESS Tail_ESS
#>         sigma_mu_1 1.030 0.120 0.837 1.287 1.001      373      138
#>         sigma_mu_2 1.296 0.144 1.063 1.659 1.001      346      231
#>    sigma_phi(1)_11 0.104 0.058 0.023 0.236 1.337        5       23
#>    sigma_phi(1)_12 0.226 0.052 0.127 0.331 1.075       32      150
#>    sigma_phi(1)_21 0.168 0.041 0.094 0.259 1.026       74      118
#>    sigma_phi(1)_22 0.192 0.043 0.114 0.279 1.017       71      100
#>  sigma_ln.sigma2_1 0.152 0.082 0.035 0.336 1.382        5       58
#>  sigma_ln.sigma2_2 0.260 0.091 0.071 0.419 1.018       45       86
#> 
#> Random Effects Correlations:
#>                            Mean    SD   2.5%  97.5%  Rhat Bulk_ESS Tail_ESS
#>                mu_1.mu_2 -0.118 0.126 -0.384  0.132 1.008      385      333
#>           mu_1.phi(1)_11 -0.073 0.284 -0.637  0.479 1.026      228      194
#>           mu_1.phi(1)_12  0.233 0.175 -0.125  0.557 1.009      245      353
#>           mu_1.phi(1)_21 -0.069 0.210 -0.472  0.333 1.006      205      269
#>           mu_1.phi(1)_22 -0.089 0.202 -0.465  0.309 1.005      322      329
#>         mu_1.ln.sigma2_1  0.105 0.284 -0.383  0.701 1.020      144      169
#>         mu_1.ln.sigma2_2  0.237 0.218 -0.222  0.634 1.030      159      221
#>           mu_2.phi(1)_11 -0.096 0.314 -0.661  0.571 1.056      129       63
#>           mu_2.phi(1)_12  0.097 0.175 -0.236  0.436 0.997      238      248
#>           mu_2.phi(1)_21  0.096 0.208 -0.302  0.500 1.012      239      342
#>           mu_2.phi(1)_22 -0.413 0.167 -0.711 -0.096 1.001      197      340
#>         mu_2.ln.sigma2_1 -0.084 0.258 -0.552  0.420 1.009      107      128
#>         mu_2.ln.sigma2_2 -0.025 0.204 -0.450  0.382 1.005      240      220
#>      phi(1)_11.phi(1)_12  0.046 0.284 -0.468  0.597 1.021      110      105
#>      phi(1)_11.phi(1)_21  0.116 0.316 -0.499  0.727 1.019       88      220
#>      phi(1)_11.phi(1)_22  0.220 0.297 -0.435  0.713 1.030       84      169
#>    phi(1)_11.ln.sigma2_1  0.217 0.324 -0.430  0.765 1.139       11      133
#>    phi(1)_11.ln.sigma2_2 -0.023 0.307 -0.567  0.541 1.000       98      259
#>      phi(1)_12.phi(1)_21 -0.296 0.231 -0.710  0.209 1.025      165      226
#>      phi(1)_12.phi(1)_22  0.159 0.231 -0.295  0.612 1.019      125      301
#>    phi(1)_12.ln.sigma2_1  0.300 0.307 -0.350  0.789 1.066       44      187
#>    phi(1)_12.ln.sigma2_2  0.053 0.274 -0.463  0.576 1.034       86      128
#>      phi(1)_21.phi(1)_22  0.139 0.253 -0.396  0.574 1.006      121      242
#>    phi(1)_21.ln.sigma2_1 -0.101 0.286 -0.632  0.463 1.013      139      240
#>    phi(1)_21.ln.sigma2_2 -0.267 0.268 -0.707  0.298 1.048      100      182
#>    phi(1)_22.ln.sigma2_1  0.278 0.283 -0.316  0.750 1.015      141      215
#>    phi(1)_22.ln.sigma2_2 -0.055 0.246 -0.558  0.394 1.012      124      326
#>  ln.sigma2_1.ln.sigma2_2 -0.114 0.292 -0.660  0.457 1.007      168      201
#> 
#> Samples were drawn using NUTS on Fri Oct 16 13:12:39 2026.
#> For each parameter, Bulk_ESS and Tail_ESS are measures of effective
#> sample size, and Rhat is the potential scale reduction factor
#> on split chains (at convergence, Rhat = 1).
# }
```
