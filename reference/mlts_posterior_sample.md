# Generate Posterior Predictive Samples for Multilevel Latent Time Series Models

The `mlts_posterior_sample()` function generates replicated datasets
from a fitted `mlts` model using draws from the posterior distribution.
The function can simulate data under the population model or based on
individual-specific (random effect) parameters.

## Usage

``` r
mlts_posterior_sample(
  fit,
  draw_person_pars = FALSE,
  n_draws = 10,
  draws = NULL,
  as_matrix = TRUE
)
```

## Arguments

- fit:

  An object of class `mlts.fit`, as returned by a fitted model using
  [`mlts()`](https://munchfab.github.io/mlts/reference/mlts-package.md).

- draw_person_pars:

  Logical. If `TRUE`, samples are generated using person-specific
  parameters (random effects). If `FALSE`, only population-level
  parameters are used. Defaults to `FALSE`.

- n_draws:

  Integer. Number of posterior draws to use for simulating replicated
  datasets. Ignored if `draws` is provided. Defaults to 10.

- draws:

  Optional integer vector indicating specific posterior draw indices to
  use. If `NULL`, `n_draws` draws are chosen with the maximum distance
  between posterior samples.

- as_matrix:

  Logical. Return replications of each variable as a matrix with
  `n_draw` rows, ready to run graphical posterior predictive checks
  using the `bayesplot` package.

## Value

A list of replicated datasets, each as a `data.frame` with columns:

- `Y_rep`:

  Replication number.

- `ID`:

  Subject/cluster ID.

- `time`:

  Time point.

- ...:

  One column per time-series variable defined in the model.

## Details

The function extracts posterior samples of population-level (and
optionally individual-level) parameters from a fitted `mlts` model and
simulates replicated datasets from the posterior predictive
distribution. Each replication corresponds to a different posterior draw
and reflects uncertainty in the model's parameters. See `PPC` for an
overview on graphical posterior predictive checks and how they can be
performed.

If `draw_person_pars = TRUE`, the function uses sampled person-specific
random effects and covariate effects from the posterior to generate new
data at the individual level. This requires that the model was fitted
with `monitor_person_pars = TRUE` in
[`mlts_fit`](https://munchfab.github.io/mlts/reference/mlts_fit.md). If
this condition is not met, the function will throw an error.

Posterior draws are either selected with the maximum distance between
posterior samples (`n_draws`) or specified manually using the `draws`
argument. Optionally, left or right censoring is respected in the
simulated data if such constraints were present in the model.

## See also

[`mlts_pp_check`](https://munchfab.github.io/mlts/reference/mlts_pp_check.md)
for plotting posterior predictive checks.

## Examples

``` r
if (FALSE) { # \dontrun{
# build a simple vector-autoregressive mlts model with two time-series variables
var_model <- mlts_model(q = 2)

# simulate data from this model with default true values
# (true values are randomly drawn from normal distribution)
var_data <- mlts_sim(
  model = var_model,
  N = 50, TP = 30, # number of units and number of measurements per unit
  default = TRUE # use default parameter values
)

# fit model
fit <- mlts_fit(
  model = var_model,
  data = var_data,
  id = "ID", ts = c("Y1", "Y2"),
  time = "time",
  monitor_person_pars = TRUE
)

# Simulate 20 replications from the posterior
yreps <- mlts_posterior_sample(fit = fit, n_draws = 20)

# Include person-specific parameters in simulation
yreps <- mlts_posterior_sample(fit = fit, draw_person_pars = TRUE)

# Use specific posterior draws
yreps <- mlts_posterior_sample(fit = fit, draws = c(10, 50, 100))
} # }
```
