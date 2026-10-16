# Changelog

## mlts (development version)

### mlts 2.0.1

#### Bug Fixes

- Selection of true values for data generation in
  [`mlts_sim()`](https://munchfab.github.io/mlts/reference/mlts_sim.md)
  was spurious for residual SDs of outcome variables. Now set to 0.5.
- Between-level standardization using
  [`mlts_standardized()`](https://munchfab.github.io/mlts/reference/mlts_standardized.md)
  failed for models with multiple time-invariant covariates.
- [`mlts_paths()`](https://munchfab.github.io/mlts/reference/mlts_paths.md)
  did not display time-invariant covariates used as outcome predictors

### mlts 2.0.0

#### New Features

- Deprecated
  [`mlts_model_paths()`](https://munchfab.github.io/mlts/reference/mlts_model_paths.md)
  and replaced with
  [`mlts_paths()`](https://munchfab.github.io/mlts/reference/mlts_paths.md).  
- Inclusion of non-lagged (i.e., same time) directed effects via the
  `incl_t0_effects` argument in
  [`mlts_model()`](https://munchfab.github.io/mlts/reference/mlts_model.md).
  Note that model formulas is not supported for these types of models at
  this point.
- Inclusion of interaction effects on the dynamic within-level via the
  `incl_interaction_effects` argument in
  [`mlts_model()`](https://munchfab.github.io/mlts/reference/mlts_model.md).
  Note that obtaining model formulas is not supported for these types of
  models at this point.
- Multiple group analyses by specifying the `group` argument in
  [`mlts_model()`](https://munchfab.github.io/mlts/reference/mlts_model.md).
- Loosen the previous restriction of allowing only bivariate VAR models
  with random innovation covariances. Models are still limited to the
  incorporation of only one latent factor capturing the innovation
  covariance, however, additional constructs can be added under the
  assumption of independent innovations.  
- Added option to include exogenous time-varying covariates (i.e.,
  without latent-mean centering) using the `is_exogenous` argument in
  [`mlts_model()`](https://munchfab.github.io/mlts/reference/mlts_model.md)
- Add option to estimate models that assume censoring on the manifest
  indicators using the `censor_left` or `censor_right` arguments in
  [`mlts_model()`](https://munchfab.github.io/mlts/reference/mlts_model.md).  
- Provide option to introduce equality constraints on loading parameters
  across levels in
  [`mlts_model()`](https://munchfab.github.io/mlts/reference/mlts_model.md)
  using `equal_loads_levels = TRUE` (i.e., fix between-level loading
  parameter of an indicator to match the respective within-level loading
  parameter).  
- Obtain posterior predictive samples using
  [`mlts_posterior_sample()`](https://munchfab.github.io/mlts/reference/mlts_posterior_sample.md)
- Run posterior predictive checks
  [`mlts_pp_check()`](https://munchfab.github.io/mlts/reference/mlts_pp_check.md)
  (still considered developmental).  
- Option to limit the number of consecutive missing values in the a
  subjects’ time series to avoid excessive numbers of `NA`s thereby
  improving model estimation times.
- [`mlts_fit()`](https://munchfab.github.io/mlts/reference/mlts_fit.md)
  now runs a small check to ensure that observed grand-means of
  `ts`-variables match prior distributions of the respective fixed
  effects. For multiple-indicator models a simple message is printed if
  any mean values fall outside the range of -10 to 10.

#### Changes

- Shortened column names in `summary` output
- Reordered sections in `summary` output
- Redundant constraint removed from re-transformed (i.e., exponentiated)
  individual variances in stan models
- When using `tinterval` in
  [`mlts_fit()`](https://munchfab.github.io/mlts/reference/mlts_fit.md)
  now slices rows with `NA`s on all `ts`-variables prior to the creation
  of the time grid
- True values of population parameters using
  [`mlts_sim()`](https://munchfab.github.io/mlts/reference/mlts_sim.md)
  with `default = TRUE` changed slightly.

### mlts 1.0.0

- Initial CRAN submission.
