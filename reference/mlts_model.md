# Build a multilevel latent time series model

Build a multilevel latent time series model

## Usage

``` r
mlts_model(
  class = c("VAR"),
  q,
  p = NULL,
  max_lag = c(1, 2, 3),
  btw_factor = TRUE,
  btw_model = NULL,
  equal_loads_levels = FALSE,
  fix_dynamics = FALSE,
  fix_inno_vars = FALSE,
  fix_inno_covs = TRUE,
  inno_covs_zero = FALSE,
  inno_covs_dir = NULL,
  fixef_zero = NULL,
  ranef_zero = NULL,
  ranef_pred = NULL,
  out_pred = NULL,
  out_pred_add_btw = NULL,
  group = NULL,
  is_exogenous = NULL,
  incl_t0_effects = NULL,
  incl_interaction_effects = NULL,
  censor_left = NULL,
  censor_right = NULL,
  silent = FALSE
)
```

## Arguments

- class:

  Character. Indicating the model type to be specified. For now
  restricted to `VAR`, the default. Future package releases might
  include additional model types.

- q:

  Integer. The number of time-varying constructs.

- p:

  Integer. For multiple-indicator models, specify a vector of length `q`
  with the number of manifest indicators per construct. If all
  constructs are measured with the same number of indicators, a single
  value is sufficient.

- max_lag:

  Integer. The maximum lag of the autoregressive effect to be included
  in the model. The maximum is 3. Defaults to 1.

- btw_factor:

  Logical. If `TRUE` (the default), a common between-level factor is
  modeled across all indicator variables per construct `q`. If `FALSE`,
  instead of a between-level factor, indicator mean levels will be
  included as individual (random) effects drawn from a joint
  multivariate normal distribution.

- btw_model:

  A list to indicate for which manifest indicator variables a common
  between-level factor should be modeled (see Details for detailed
  instructions). At this point restricted to one factor per latent
  construct.

- equal_loads_levels:

  Logical. For multiple-indicator model with `btw_factor = TRUE`, if
  `TRUE`, factor loadings of the same indicators are assumed to be equal
  across levels. Note, that the first indicator loading parameters
  remain fixed to `1`.

- fix_dynamics:

  Logical. Fix all random effect variances of autoregressive and
  cross-lagged effects to zero (constraining parameters to be equal
  across clusters).

- fix_inno_vars:

  Logical. Fix all random effect variances of innovation variances to
  zero (constraining parameters to be equal across clusters).

- fix_inno_covs:

  Logical. Fix all random effect variances of innovation covariances to
  zero (constraining parameters to be equal across clusters).

- inno_covs_zero:

  Logical. Set to `TRUE` to treat all innovations as independent.

- inno_covs_dir:

  For bivariate VAR models with person-specific innovation covariances,
  a latent variable approach is applied (for a detailed description, see
  Hamaker et al., 2018). by specifying an additional factor that loads
  onto the contemporaneous innovations of both constructs, capturing the
  shared variance of innovations, that is not predicted by the previous
  time points. The loading parameters of this latent factor, however,
  have to be restricted in accordance with researchers assumptions about
  the sign of the association between innovations across construct.
  Hence, if innovations at time \$t\$ are assumed to be positively
  correlated across clusters, set the argument to `pos`, or `neg`
  respectively.

- fixef_zero:

  Character. A character vector to index which fixed effects (referring
  to the parameter labels in `model$Param`) should be constrained to
  zero (Note: this also results in removing the random effect variance
  of the respective parameter).

- ranef_zero:

  Character. A character vector to index which random effect variances
  (referring to the parameter labels in `model$Param`) should be
  constrained to zero.

- ranef_pred:

  A character vector or a named list. Include between-level covariate(s)
  as predictor(s) of all random effects in `model` by entering a vector
  of unique variable names. Alternatively, to include between-level
  covariates or differing sets of between-level covariates as predictors
  of specific random effects, a named list (using the labels in
  `model$Param`) can be entered (see examples). Note that if a named
  list is provided, all names that do not match random parameters in
  `model` will be ignored. Note that variables entered in `ranef_pred`
  will be grand-mean centered by default when fitting the model with
  `mlts_fit`.

- out_pred:

  A character vector or a named list. Include between-level outcome(s)
  to be regressed on all random effects in `model` by entering a vector
  of unique variable names. Alternatively, to include multiple
  between-level outcomes regressed differing sets of specific random
  effects, a named list (using the labels in `model$Param`) can be
  entered (see examples). Note that if a named list is provided, all
  character strings in the vector of each list (with independent
  variables) element that do not match random effect parameter names in
  `model$Param` will be treated as additional between-level predictors.

- out_pred_add_btw:

  A character vector. If `out_pred` is a character (vector), all inputs
  will be treated as between-level covariates to be used as additional
  predictors of all outcomes specified in `out_pred`.

- group:

  An integer specifying the number of groups (not yet supported). Add a
  binary coded (0 vs. 1) variable to include group differences in fixed
  effects (intercepts). When dynamic or variance parameters are allowed
  to vary by cluster, you can enter the grouping variable to `re_pred`.

- is_exogenous:

  Integer or a vector of integers. Indicate if any of the constructs
  should be treated as exogenous (i.e., no latent mean centering will be
  performed). Probable use case: Adding a dichotomous time-varying
  predictor variable.

- incl_t0_effects:

  A character vector. Experimental: Add contemporaneous effects to the
  model. For example, to include an effect of the first construct on the
  second construct at time \$t\$, following the general pattern for
  naming of dynamic parameters in the mlts framework, can be included by
  specifying `phi(0)_21` where the `0` indicates the lag, the first
  subscript letter (`2`) the dependent, and the latter subscript (`1`)
  the independent construct. The respective within-level
  correlation/covariance of innovations between involved constructs will
  be excluded from the model accordingly.

- incl_interaction_effects:

  A character vector. Add interaction terms on the dynamic within-level.
  For example, to add an interaction term between the first construct at
  time \$t\$ (lag of 0) and the second construct at \$t-1\$ (lag of 1)
  to the prediction of the second construct at time \$t\$ specify
  `incl_interaction_effects = phi(i)_2.2(1)1(0)`. where the `i`
  indicates an interaction effect, the first subscript letter (`2`) the
  dependent, and the latter subscripts after the dot (i.e., `2(1)` and
  `1(0)`) the independent constructs involved in the interaction each
  followed by the respective lag in brackets. Note, that in this case
  the respective lag 0 effects need to be included separately using
  `incl_t0_effects`.

- censor_left:

  Numeric. If an input is provided (i.e., a single numeric value) a
  left-censored version of the model will be estimated by treating all
  observations (of manifest indicators) at the censoring threshold
  (i.e., usually the lower bound of the scale) to be treated as missing
  during model estimation. These missing values (observations at the
  value of `censor_left`) are replaced with imputed values (declared as
  parameters in the stan model) with an upper limit of `censor_left`
  (see
  https://mc-stan.org/docs/stan-users-guide/truncation-censoring.html).
  Note that all manifest variables are affected by the censoring. To
  prevent individual variables from being treated as censored you could
  change the scale of the respective variable(s) so that all values
  exceed the censoring threshold.

- censor_right:

  Numeric. Developmental. Similar to `censor_left` but assumes variables
  to be censored on the upper bound of the scale. Can be combined with
  `censor_left`.

- silent:

  logical. Set to `TRUE` to suppress warnings and messages.

## Value

An object of class `data.frame` with the following columns:

- Model:

  Indicates if the parameter in the respective row is part of the
  structural, or the measurement model (if multiple indicators per
  construct are provided)

- Level:

  Parameter on the between- or within-level.

- Type:

  Describes the parameter type.

- Param:

  Parameter names to be referred to in arguments of `mlts_model`.

- Param_Label:

  Parameter labels (additional option to address specific parameters).

- isRandom:

  Indicates which within-level parameters are modeled as random (1) or a
  constant across clusters (0).

- Constraint:

  Optional. Included if multiple-indicators per construct (p \> 1) are
  provided. Constraints on measurement model parameters can be changed
  by overwriting the respective value in `model`. Possible inputs are
  "free", "= 0" (for SDs of measurement error variances), and "= 1" (for
  loading parameters).

- prior_type:

  Contains the parameters' prior distribution used in `mlts_fit` (prior
  classes can not be changed at this point).

- prior_location:

  Location values of the parameters' prior distribution used in
  `mlts_fit` (can be changed to any real value by overwriting the
  respective value in `model`).

- prior_scale:

  Scale values of the parameters' prior distribution used in `mlts_fit`
  (can be changed to any real value by overwriting the respective value
  in `model`).

## References

Hamaker, E. L., Asparouhov, T., Brose, A., Schmiedek, F., & Muthén, B.
(2018). At the frontiers of modeling intensive longitudinal data:
Dynamic structural equation models for the affective measurements from
the COGITO study. *Multivariate behavioral research*, *53*(6), 820-841.
[doi:10.1080/00273171.2018.1446819](https://doi.org/10.1080/00273171.2018.1446819)

## Examples

``` r
# \donttest{
 # To illustrate the general model building procedure, starting with a simple
 # two-level AR(1) model with person-specific individual means, AR effects,
 # and innovation variances (the default option when using mlts_model() and q = 1).
 model <- mlts_model(q = 1)

 # All model parameters (with their labels stored in model$Param) can be inspected by calling:
 model
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

 # Possible model extensions/restrictions:
 # 1. Introducing additional parameter constraints, such as fixing specific
 #    parameters to a constant value by setting the respective random effect
 #    variances to zero, such as e.g. (log) innovation variances
 model <- mlts_model(q = 1, ranef_zero = "ln.sigma2_1")
 #    Note that setting the argument `fix_inno_vars` to `TRUE` provides
 #    a shortcut to fixing the innovation variances of all constructs
 #    (if q >= 1) to a constant.

 # 2. Including a multiple indicator model, where the construct is measured by
 #    multiple indicators (here, p = 3 indicators)
 model <- mlts_model(
          q = 1, # the number of time-varying constructs
          p = 3, # the number of manifest indicators
          # assuming a common between-level factor (the default)
          btw_factor = TRUE
        )

 # 3. Incorporating between-level variables. For example, inclusion of
 #    an additional between-level variable ("cov1") as predictor of all
 #    (ranef_pred = "cov1") or a specific set of random effects
 #    (ranef_pred = list("phi(1)_11") = "cov1"), an external outcome (e.g., "out1")
 #    to be predicted by all (out_pred = "out1") or specific random effects
 #    (out_pred = list("out1" = c("etaB_1", "phi(1)_11")), using the latent
 #    between-level factor trait scores (etaB_1) and individual first-order
 #    autoregressive effects (phi(1)_11) as joint predictors of outcome "out1".
 model <- mlts_model(
            q = 1,
            p = 3,
            fix_inno_vars = TRUE,
            ranef_pred = "cov1",
            out_pred = list("out1" = c("etaB_1", "phi(1)_11"))
           )
 #    Note that the names of the random effect parameters must match the
 #    parameter labels provided in model$Param, the result of the
 #    mlts_model()-functions.

# }
```
