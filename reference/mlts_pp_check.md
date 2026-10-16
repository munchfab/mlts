# Posterior Predictive Checks for Multilevel Latent Time Series Models

Developemental This function plots posterior predictive distributions of
one or multiple fitted `mlts.fit` models. Simulated data from posterior
draws are compared to the observed data to visually assess model fit.

## Usage

``` r
mlts_pp_check(
  fit,
  fit_list = NULL,
  ts = NULL,
  y_reps = NULL,
  by_cluster = FALSE,
  by_group = FALSE,
  cluster_ids = NULL,
  draw_person_pars = FALSE,
  n_draws = 10,
  draws = NULL,
  add_y_obs = TRUE,
  model_lab = NULL,
  y_rep_col = NULL,
  y_obs_col = "#009E73",
  y_obs_lw = 1.1,
  y_rep_lw = 0.5,
  y_rep_alpha = 0.5
)
```

## Arguments

- fit:

  A fitted model object of class `mlts.fit`. Only used if `fit_list` is
  `NULL`.

- fit_list:

  An optional list of fitted `mlts.fit` objects for model comparison. If
  provided, `fit` is ignored.

- ts:

  Optional vector of variable names to include in the plot.

- y_reps:

  Optional. A list of posterior predictive samples (as returned by
  [`mlts_posterior_sample`](https://munchfab.github.io/mlts/reference/mlts_posterior_sample.md))
  for a single model. If `NULL`, samples are generated within the
  function.

- by_cluster:

  Logical. If `TRUE`, density plots are faceted by individual and
  time-series variable. If `FALSE`, only time-series variables are used
  for faceting. Default is `FALSE`.

- by_group:

  Logical. If `TRUE`, density plots are faceted by grouping and
  time-series variable. If `FALSE`, only time-series variables are used
  for faceting. Default is `FALSE`.

- cluster_ids:

  Optional vector of cluster IDs to include in the plot. If `NULL`, all
  IDs are shown.

- draw_person_pars:

  Logical. If `TRUE`, samples are generated using person-specific
  parameters (random effects). If `FALSE`, only population-level
  parameters are used. Defaults to `FALSE`.

- n_draws:

  Integer. Number of posterior draws to use for generating replicated
  datasets. Defaults to 20. Ignored if `draws` is specified.

- draws:

  Optional vector of indices specifying which posterior draws to use. If
  `NULL`, `n_draws` samples are drawn randomly.

- add_y_obs:

  Logical. Whether to include the observed data distribution in the
  plot. Defaults to `TRUE`.

- model_lab:

  Optional character vector with labels for each model in `fit_list`. If
  `NULL`, defaults to "Model 1", "Model 2", etc.

- y_rep_col:

  Optional vector of colors for the posterior predictive densities of
  each model. If `NULL`, a default color palette is used.

- y_obs_col:

  Color for the observed data distribution. Default is `"#009E73"`.

- y_obs_lw:

  Line width of the observed data density curve. Default is `1.1`.

- y_rep_lw:

  Line width of the observed data density curve. Default is `0.5`.

- y_rep_alpha:

  Alpha transparency for the predictive density curves. Default is
  `0.5`.

## Value

A `ggplot` object showing density curves of observed and replicated data
across time-series variables (and optionally across individuals).

## Details

This function performs graphical posterior predictive checks by
overlaying kernel density estimates of replicated data from the
posterior with the observed data. This can be used to visually assess
how well a fitted model captures key distributional aspects of the
observed time series. If `fit_list` is specified, multiple models can be
compared side-by-side in the same plot.

If `draw_person_pars = TRUE`, simulated datasets incorporate
subject-specific effects (random effects). This requires that
`monitor_person_pars = TRUE` was set during model fitting.

## See also

[`mlts_posterior_sample`](https://munchfab.github.io/mlts/reference/mlts_posterior_sample.md)
for generating replicated data samples.

## Examples

``` r
if (FALSE) { # \dontrun{
# Set up AR(1) model
ar1 <- mlts_model(q = 1, censor_left = -1)

# Simulate data under the AR(1) model
simData <- mlts_sim(model = ar1, N = 50, TP =100, default = TRUE)

# Fit the model
fit_AR <- mlts_fit(model = ar1, data = simData$data,
                   id = "ID", ts = "Y1", monitor_person_pars = TRUE)

# Run posterior predictive check
mlts_pp_check(fit = fit_AR,
              model_lab = "AR(1)",
              y_rep_col = "steelblue")

} # }
```
