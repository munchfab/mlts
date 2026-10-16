# Include Predictors for Random Effects on the Between Level

``` r

library(mlts)
```

For this example, we use data from the website
[PhysioNet](https://physionet.org/) (Goldberger et al., 2000), which
hosts open data for physiological research. The data come from a study
by Hausdorff et al. (1996), who examined walking stride intervals (gait)
of 15 subjects (5 healthy young adults, age 23 – 29; 5 healthy old
adults, age 71 – 77; and 5 older adults with Parkinson’s disease). On
the [data website](https://physionet.org/content/gaitdb/1.0.0/), the
description regarding walking stride intervals read:

> Subjects walked continuously on level ground around an obstacle-free
> path. The stride interval was measured using ultra-thin, force
> sensitive resistors placed inside the shoe. The analog force signal
> was sampled at 300 Hz with a 12 bit A/D converter, using an
> ambulatory, ankle-worn microcomputer that also recorded the data.
> Subsequently, the time between foot-strikes was automatically
> computed. The method for determining the stride interval is a
> modification of a previously validated method that has been shown to
> agree with force-platform measures, a “gold” standard. Data were
> collected from the healthy subjects as they walked in a roughly
> circular path for 15 minutes, and from the subjects with Parkinson’s
> disease as they walked for 6 minutes up and down a long hallway.

We will model subjects’ walking stride intervals in an autoregressive
model. Furthermore, we will examine if subjects’ health status (healthy
or Parkinson’s disease) can explain variation in random model parameters
(i.e., mean stride interval $`\mu_1`$, autoregressive effect
$`\phi_{(1)11}`$, or log innovation variance
$`\ln(\sigma^2_{\zeta_{1}})`$).

``` r

load("gait_data.rda")
head(gait_data)
#>   subject age       group   time stride_interval healthy
#> 1       1  76 healthy_old 30.797           1.023       1
#> 2       1  76 healthy_old 31.820           1.030       1
#> 3       1  76 healthy_old 32.850           1.017       1
#> 4       1  76 healthy_old 33.867           1.027       1
#> 5       1  76 healthy_old 34.893           1.043       1
#> 6       1  76 healthy_old 35.937           1.027       1
```

The variables included in this data set are

- `subject`: the unit identifier (ID)
- `age`: subjects’ age in years
- `group`: group variable with levels `healthy_old`, `healthy_young`,
  and `pd_old` (for old adults with Parkinson’s disease)
- `time`: the variable containing the measurement point after begin of
  the study, in seconds. It starts after a brief “warm-up” period
- `stride_interval`: subjects’ stride interval in seconds (i.e., the
  time between heel strikes)

## Autoregressive model without predictors

First, we will fit a simple autoregressive model for the stride interval
time series. The model setup is the same as in the vignette “A Simple
Example: Multilevel Manifest AR(1) Model”. We set it up with

``` r

ar1_model1 <- mlts_model(q = 1)
```

We can check the parameters present in the model by just calling the
object:

``` r

ar1_model1
#>                              Model   Level             Type
#> mu_1                    Structural  Within     Fixed effect
#> phi(1)_11               Structural  Within     Fixed effect
#> ln.sigma2_1             Structural  Within     Fixed effect
#> sigma_mu_1              Structural Between Random effect SD
#> sigma_phi(1)_11         Structural Between Random effect SD
#> sigma_ln.sigma2_1       Structural Between Random effect SD
#> r_mu_1.phi(1)_11        Structural Between   RE correlation
#> r_mu_1.ln.sigma2_1      Structural Between   RE correlation
#> r_phi(1)_11.ln.sigma2_1 Structural Between   RE correlation
#>                                           Param             Param_Label
#> mu_1                                       mu_1                   Trait
#> phi(1)_11                             phi(1)_11                 Dynamic
#> ln.sigma2_1                         ln.sigma2_1 Log Innovation Variance
#> sigma_mu_1                           sigma_mu_1                   Trait
#> sigma_phi(1)_11                 sigma_phi(1)_11                 Dynamic
#> sigma_ln.sigma2_1             sigma_ln.sigma2_1 Log Innovation Variance
#> r_mu_1.phi(1)_11               r_mu_1.phi(1)_11                  RE Cor
#> r_mu_1.ln.sigma2_1           r_mu_1.ln.sigma2_1                  RE Cor
#> r_phi(1)_11.ln.sigma2_1 r_phi(1)_11.ln.sigma2_1                  RE Cor
#>                         isRandom prior_type prior_location prior_scale
#> mu_1                           1     normal              0        10.0
#> phi(1)_11                      1     normal              0         2.0
#> ln.sigma2_1                    1     normal              0        10.0
#> sigma_mu_1                     0     cauchy              0         2.5
#> sigma_phi(1)_11                0     cauchy              0         2.5
#> sigma_ln.sigma2_1              0     cauchy              0         2.5
#> r_mu_1.phi(1)_11               0        LKJ              1          NA
#> r_mu_1.ln.sigma2_1             0        LKJ              1          NA
#> r_phi(1)_11.ln.sigma2_1        0        LKJ              1          NA
```

For convenience, we can furthermore check the associated path model:

``` r

mlts_model_paths(ar1_model1)
```

![](../articles/pathmodel_ar1.png)

We fit the model with
[`mlts_fit()`](https://munchfab.github.io/mlts/reference/mlts_fit.md).
We need to provide the unit identifier and the time-series construct we
wish to examine:

``` r

ar1_fit1 <- mlts_fit(
  model = ar1_model1,
  data = gait_data,
  id = "subject",
  ts = "stride_interval",
  monitor_person_pars = TRUE,
  iter = 1000
)
```

``` r

summary(ar1_fit1)
```

    #> Call:
    #> mlts_model(q = 1, max_lag = 1)
    #> Time series variables as indicated by parameter subscripts: 
    #>    1 --> stride_interval
    #> Data: 9144 observations in 15 IDs
    #> Model convergence criteria: 
    #>   Maximum Potential Scale Reduction Factor (PSR; Rhat): 1.006 (should be < 1.01)
    #>   Minimum Bulk ESS: 482 (should be > 200, 100 per chain) 
    #>   Minimum Tail ESS: 506 (should be > 200, 100 per chain) 
    #>   Number of divergent transitions: 0 (should be 0) 
    #> 
    #> Fixed Effects:
    #>              Post. Mean Post. Median Post. SD   2.5%  97.5%  Rhat Bulk_ESS
    #>         mu_1      1.096        1.096    0.033  1.029  1.157 1.004      747
    #>    phi(1)_11      0.348        0.348    0.080  0.185  0.512 1.006      723
    #>  ln.sigma2_1     -6.958       -6.959    0.397 -7.777 -6.216 1.003      622
    #>  Tail_ESS
    #>       713
    #>       506
    #>       587
    #> 
    #> Random Effects SDs:
    #>              Post. Mean Post. Median Post. SD  2.5% 97.5%  Rhat Bulk_ESS
    #>         mu_1      0.132        0.128    0.028 0.091 0.200 1.001      778
    #>    phi(1)_11      0.296        0.287    0.065 0.200 0.450 1.003      664
    #>  ln.sigma2_1      1.553        1.509    0.310 1.073 2.275 1.002      678
    #>  Tail_ESS
    #>       785
    #>       635
    #>       681
    #> 
    #> Random Effects Correlations:
    #>                        Post. Mean Post. Median Post. SD   2.5%  97.5%  Rhat
    #>         mu_1.phi(1)_11      0.015        0.018    0.246 -0.463  0.491 1.002
    #>       mu_1.ln.sigma2_1      0.480        0.499    0.197  0.058  0.801 1.001
    #>  phi(1)_11.ln.sigma2_1     -0.459       -0.487    0.196 -0.774 -0.029 1.002
    #>  Bulk_ESS Tail_ESS
    #>       714      685
    #>       554      695
    #>       539      718
    #> 
    #> Samples were drawn using NUTS on Thu May 16 15:17:38 2024.
    #> For each parameter, Bulk_ESS and Tail_ESS are measures of effective
    #> sample size, and Rhat is the potential scale reduction factor
    #> on split chains (at convergence, Rhat = 1).

The fixed effect of the subject mean of stride interval is estimated at
1.096. Individual mean stride intervals fluctuate around this fixed
effect with a standard deviation of 0.132 (see `Random Effect SDs`). The
fixed effect of the first-order autoregressive effect (see `phi(1)_11`
in the section `Fixed Effects`) is estimated at 0.348, with random
effects standard deviation of 0.296. The fixed effect of log innovation
variance is estimated at -6.958. To be on the original scale, we have to
exponentiate it: exp(0.132) = 0.001. The random effects’ standard
deviation of log innovation variance is estimated at 1.553.

There is also a substantial correlation between random effects of the
person mean (`mu_1`) and log innovation variance (`ln.sigma2_1`) of
0.48, indicating that subjects with a longer stride interval also
display higher variation in their stride interval. Furthermore, there is
a substantial negative correlation of -0.459 between the autoregressive
effect `phi(1)_11` and log innovation variance, indicating that subjects
displaying a higher stability in their stride interval tend to display
lower variation and vice versa.

## Include a predictor for random effects

We can now try to explain differences in model parameters by another
covariate. In this example, we will use the subjects health status as
explanatory variable, which is binary in our case (i.e., `1` for healthy
subjects and `0` for subjects with Parkinson’s disease). We will include
this covariate on the *between-level* to explain random effect
variation, thus it has to be stable *within* subjects. Note furthermore
that the variable has to be numeric or integer for Stan to be
acceptable.

We set this model up with `mtls_model()`. To predict all random effects
present in the model, we can just provide the variable name in the
argument `ranef_pred`:

``` r

ar1_model2 <- mlts_model(q = 1, ranef_pred = "healthy")
```

We can see that three new parameters are added on the between level,
specifying the regression weights for predicition of random effects
(i.e., `b_mu_1.ON.healthy` is the regression weight of the predictor
`healthy` on the dependent variable `mu_1`):

``` r

ar1_model2
#>                               Model   Level             Type
#> mu_1                     Structural  Within     Fixed effect
#> phi(1)_11                Structural  Within     Fixed effect
#> ln.sigma2_1              Structural  Within     Fixed effect
#> sigma_mu_1               Structural Between Random effect SD
#> sigma_phi(1)_11          Structural Between Random effect SD
#> sigma_ln.sigma2_1        Structural Between Random effect SD
#> r_mu_1.phi(1)_11         Structural Between   RE correlation
#> r_mu_1.ln.sigma2_1       Structural Between   RE correlation
#> r_phi(1)_11.ln.sigma2_1  Structural Between   RE correlation
#> b_mu_1.ON.healthy        Structural Between    RE prediction
#> b_phi(1)_11.ON.healthy   Structural Between    RE prediction
#> b_ln.sigma2_1.ON.healthy Structural Between    RE prediction
#>                                             Param             Param_Label
#> mu_1                                         mu_1                   Trait
#> phi(1)_11                               phi(1)_11                 Dynamic
#> ln.sigma2_1                           ln.sigma2_1 Log Innovation Variance
#> sigma_mu_1                             sigma_mu_1                   Trait
#> sigma_phi(1)_11                   sigma_phi(1)_11                 Dynamic
#> sigma_ln.sigma2_1               sigma_ln.sigma2_1 Log Innovation Variance
#> r_mu_1.phi(1)_11                 r_mu_1.phi(1)_11                  RE Cor
#> r_mu_1.ln.sigma2_1             r_mu_1.ln.sigma2_1                  RE Cor
#> r_phi(1)_11.ln.sigma2_1   r_phi(1)_11.ln.sigma2_1                  RE Cor
#> b_mu_1.ON.healthy               b_mu_1.ON.healthy       regression weight
#> b_phi(1)_11.ON.healthy     b_phi(1)_11.ON.healthy       regression weight
#> b_ln.sigma2_1.ON.healthy b_ln.sigma2_1.ON.healthy       regression weight
#>                          isRandom prior_type prior_location prior_scale
#> mu_1                            1     normal              0        10.0
#> phi(1)_11                       1     normal              0         2.0
#> ln.sigma2_1                     1     normal              0        10.0
#> sigma_mu_1                      0     cauchy              0         2.5
#> sigma_phi(1)_11                 0     cauchy              0         2.5
#> sigma_ln.sigma2_1               0     cauchy              0         2.5
#> r_mu_1.phi(1)_11                0        LKJ              1          NA
#> r_mu_1.ln.sigma2_1              0        LKJ              1          NA
#> r_phi(1)_11.ln.sigma2_1         0        LKJ              1          NA
#> b_mu_1.ON.healthy               0     normal              0        10.0
#> b_phi(1)_11.ON.healthy          0     normal              0        10.0
#> b_ln.sigma2_1.ON.healthy        0     normal              0        10.0
```

For convenience, we can plot the path model (note that only the
between-level model is shown here, as decomposition and within-level
model are the same as before):

``` r

mlts_model_paths(ar1_model2)
```

![](../articles/pathmodel_betw_preds.png)

We fit the model with
[`mlts_fit()`](https://munchfab.github.io/mlts/reference/mlts_fit.md).
We need to provide the unit identifier and the time-series construct we
wish to examine. If the explanatory variable is called the same as
provided in the call to
[`mlts_model()`](https://munchfab.github.io/mlts/reference/mlts_model.md),
it does not need to be specified again. Explanatory variables for random
effects are grand-mean-centered per default. As our group variable is
dichotomous, we set the argument `center_covs` to `FALSE` so that we can
interpret the intercept of the random model parameters as mean parameter
in the group of subjects with Parkinson’s disease and the effect of the
explanatory variable as group difference. To obtain standardized
estimates for parameters on the within- and between-level, we set
`monitor_person_pars = TRUE`. For large models, this may greatly
increase sampling times, so this argument is set to `FALSE` by default.

``` r

ar1_fit2 <- mlts_fit(
  model = ar1_model2,
  data = gait_data,
  id = "subject",
  ts = "stride_interval",
  center_covs = FALSE,
  monitor_person_pars = TRUE,
  iter = 1000
)
```

``` r

summary(ar1_fit2)
```

    #> Call:
    #> mlts_model(q = 1, max_lag = 1)
    #> Time series variables as indicated by parameter subscripts: 
    #>    1 --> stride_interval
    #> Data: 9144 observations in 15 IDs
    #> Model convergence criteria: 
    #>   Maximum Potential Scale Reduction Factor (PSR; Rhat): 1.017 (should be < 1.01)
    #>   Minimum Bulk ESS: 408 (should be > 200, 100 per chain) 
    #>   Minimum Tail ESS: 362 (should be > 200, 100 per chain) 
    #>   Number of divergent transitions: 0 (should be 0) 
    #> 
    #> Fixed Effects:
    #>              Post. Mean Post. Median Post. SD   2.5%  97.5%  Rhat Bulk_ESS
    #>         mu_1      1.161        1.158    0.060  1.048  1.291 1.003      429
    #>    phi(1)_11      0.160        0.166    0.126 -0.097  0.407 1.001      441
    #>  ln.sigma2_1     -5.268       -5.263    0.391 -6.106 -4.463 1.004      460
    #>  Tail_ESS
    #>       509
    #>       487
    #>       450
    #> 
    #> Random Effects SDs:
    #>              Post. Mean Post. Median Post. SD  2.5% 97.5%  Rhat Bulk_ESS Tail_ESS
    #>         mu_1      0.129        0.124    0.029 0.088 0.202 0.999      845      743
    #>    phi(1)_11      0.267        0.257    0.064 0.171 0.428 1.003      922      663
    #>  ln.sigma2_1      0.873        0.842    0.189 0.590 1.310 1.000      862      622
    #>  Tail_ESS
    #>       743
    #>       663
    #>       622
    #> 
    #> Random Effects Correlations:
    #>                        Post. Mean Post. Median Post. SD   2.5% 97.5%  Rhat Bulk_ESS
    #>         mu_1.phi(1)_11      0.181        0.190    0.248 -0.339 0.625 1.001      710
    #>       mu_1.ln.sigma2_1      0.420        0.447    0.223 -0.054 0.781 1.000      842
    #>  phi(1)_11.ln.sigma2_1     -0.190       -0.200    0.249 -0.643 0.312 1.004      810
    #>  Tail_ESS
    #>       811
    #>       848
    #>       749
    #> 
    #> Random Effects Regressed On:
    #>                        Post. Mean Post. Median Post. SD   2.5%  97.5%  Rhat
    #>         mu_1 ~ healthy     -0.096       -0.093    0.073 -0.243  0.046 1.002
    #>    phi(1)_11 ~ healthy      0.285        0.277    0.155 -0.032  0.597 0.999
    #>  ln.sigma2_1 ~ healthy     -2.540       -2.533    0.479 -3.532 -1.580 1.003
    #>  Bulk_ESS Tail_ESS
    #>       417      423
    #>       468      452
    #>       429      362
    #> 
    #> Samples were drawn using NUTS on Fri Jun 21 11:04:08 2024.
    #> For each parameter, Bulk_ESS and Tail_ESS are measures of effective
    #> sample size, and Rhat is the potential scale reduction factor
    #> on split chains (at convergence, Rhat = 1).

The [`summary()`](https://rdrr.io/r/base/summary.html) now shows a new
section called `Random Effects Regressed On`, which shows the regression
weights of the explanatory variable `healthy`. We can see that subjects’
health status has effects on all between-level model parameters the
autoregressive effect `phi(1)_11` and log innovation variance
`ln.sigma2_1`. However, only for the effect on log innovation variance,
the 95% credible intervals of the regression weight does not include 0.
In particular, healthy subjects’ (with a value of `1` in the variable
`healthy`) have, on average, a stride interval that is -0.096 seconds
shorter than subjects with Parkinson’s disease. Furthermore, healthy
subjects display, on average, an autoregressive effect that is 0.285
units higher than for subjects with Parkinson’s disease. At last, the
log innovation variance is on average -2.54 units lower for healthy
subjects in comparison to subjects with Parkinson’s disesase.

To obtain standardized estimates of the parameters, we call
[`mlts_standardized()`](https://munchfab.github.io/mlts/reference/mlts_standardized.md)
on the `mltsfit`-object. The argument `what = "both"` controls that
standardized estimates are computed for the within- and between level
(by default, it is done for the between-level only).

``` r

mlts_standardized(ar1_fit2, what = "both")
```

    #> $`Between-level standardized`
    #>             Type                    Param Std.Est    SD   2.5%  97.5%
    #> 10 RE prediction        b_mu_1.ON.healthy  -0.326 0.220 -0.690  0.147
    #> 11 RE prediction   b_phi(1)_11.ON.healthy   0.446 0.207 -0.056  0.745
    #> 12 RE prediction b_ln.sigma2_1.ON.healthy  -0.808 0.085 -0.916 -0.580
    #> 
    #> $`Within-level standardized effects averaged over clusters`
    #>           Type     Param Std.Est    SD  2.5% 97.5%
    #> 1 Fixed effect phi(1)_11    0.35 0.011 0.329 0.373

## References

Goldberger, A. L., Amaral, L. A. N., Glass, L., Hausdorff, J. M.,
Ivanov, P. Ch., Mark, R. G., Mietus, J. E., Moody, G. B., Peng, C.-K., &
Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet:
Components of a new research resource for complex physiologic signals.
*Circulation*, *101*(23), e215–e220.
<https://doi.org/10.1161/01.CIR.101.23.e215>

Hausdorff, J., Purdon, P., Peng, C., Ladin, Z., Wei, J., & Goldberger,
A. (1996). Fractal dynamics of human gait: Stability of long-range
correlations in stride interval fluctuations. *Journal of Applied
Physiology*, *80*(5), 1448—1457.
<https://doi.org/10.1152/jappl.1996.80.5.1448>
