# Reliability estimates for individual model parameters

Computes reliability estimates for individual-specific parameters from a
fitted `mlts` model using approaches recommended for intensive
longitudinal data.

## Usage

``` r
mlts_rel(mltsfit, method = c("EAP", "RMU"), prob = 0.95, seed = NULL)
```

## Arguments

- mltsfit:

  An object returned by
  [`mlts()`](https://munchfab.github.io/mlts/reference/mlts-package.md).

- method:

  Character string specifying the reliability estimation method. Must be
  one of `"EAP"` or `"RMU"`.

- prob:

  Numeric value between 0 and 1 specifying the probability mass covered
  by the equal-tailed credible interval (ETI) in the `"RMU"` method.
  Defaults to `0.95`.

- seed:

  Integer used in the `"RMU"` method.

## Value

If `method = "EAP"`:

- A named numeric vector containing reliability estimates for all random
  individual parameters.

If `method = "RMU"`:

- A data.frame with one row per random parameter containing:

  - `rel_mean`: Mean reliability estimate

  - `rel_mdn`: Median reliability estimate

  - `rel_sd`: Standard deviation of reliability estimates

  - `rel_ci.lb`: Lower bound of the equal-tailed credible interval

  - `rel_ci.ub`: Upper bound of the equal-tailed credible interval

## Details

Reliability is conceptualized as the proportion of true individual
differences relative to total variance, accounting for estimation
uncertainty.

Two methods are currently implemented:

- `"EAP"`: Reliability based on posterior mean estimates and posterior
  standard deviations.

- `"RMU"`: Reliability based on repeated random split-half correlations
  of posterior draws, summarized by location and uncertainty.

The `"EAP"` method estimates reliability as the ratio of true score
variance to total variance, where total variance is decomposed into
between-person variability in the posterior mean estimates and the
expected posterior uncertainty (measurement error). Specifically,
reliability is defined as \$\$ Rel =
\frac{Var(\hat{\theta})}{Var(\hat{\theta}) + E(SE^2)} \$\$ where
\\\hat{\theta}\\ denotes the posterior mean (EAP) estimates of the
individual parameters and \\SE\\ their posterior standard deviations.
This formulation follows the conceptualization of reliability as the
proportion of variance in observed estimates that reflects stable
between-person differences, rather than estimation uncertainty.

The `"RMU"` method implements a resampling-based approach to quantify
the reliability of individual differences in model parameters. Posterior
draws are repeatedly split into independent subsets, and correlations
between these subsets are computed across individuals. These split-half
correlations provide an estimate of the stability (replicability) of
individual differences under posterior uncertainty. The resulting
distribution of correlations is summarized to obtain point estimates and
uncertainty intervals (equal-tailed credible intervals) for reliability.

Both approaches explicitly account for uncertainty in individual
parameter estimation and are designed for intensive longitudinal data,
where reliability pertains to the stability of individual-specific model
parameters rather than observed variables.

The `"RMU"` implementation is described in: Bignardi, G., Kievit, R. A.,
& Bürkner, P.C. (2025), *A general method for estimating reliability
using Bayesian Measurement Uncertainty*, PsyArXiv.
<https://osf.io/preprints/psyarxiv/h54k8_v1>

## Examples

``` r
if (FALSE) { # \dontrun{

# Computation of individual parameter reliability for an AR(1) model:
## Model specification
AR1 <- mlts_model(q = 1)

## Simulate data
simData <- mlts_sim(model = AR1, N = 70, TP = 70, default = TRUE)

## Estimate the model and store MCMC draws of individual parameters
fit <- mlts_fit(
  model = AR1, data = simData$data,
  id = "ID", ts = "Y1",
  monitor_person_pars = TRUE # !important
)

# Person separation reliability
mlts_rel(fit, method = "EAP")

# Relative Measurement Uncertainty (RMU)
mlts_rel(fit, method = "RMU")
} # }
```
