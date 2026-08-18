MLTS
================

<!-- badges: start -->

[![](https://cranlogs.r-pkg.org/badges/mlts)](https://cranlogs.r-pkg.org/badges/mlts)
[![CRAN
status](https://www.r-pkg.org/badges/version/mlts)](https://CRAN.R-project.org/package=mlts)

<!-- badges: end -->

The `mlts` package allows fitting multilevel manifest or latent time
series models and dynamic structural equation models (DSEM), as
described in Asparouhov et al. (2018). It relies on
[Stan](https://mc-stan.org) (Stan Development Team, 2023b) and the
`rstan` package (Stan Development Team, 2023a) for Bayesian inference.
The package is designed for researchers working with intensive
longitudinal data, e.g., data from ambulatory assessment studies or
experience sampling methods. Models allow for missing data, unequal
numbers of observations across units, unequally spaced observations, and
approximation of continuous time processes. Parameters can be subjected
to restrictions (i.e., fixations or equality constraints) and latent
variables can be incorporated using a measurement model. For the models
hypothesized by the user, a LaTex formula and path models can be
inspected, and model parameters can be plotted.

At this stage, the package incorporates the means for:

- Etimation of two-level vector autoregressive models, which flexibly
  handle

  - multiple time-varying constructs (maximum of nine),

  - lagged effects of higher order (up to lag three),

  - user-specified within-level dynamics (unique lagged effects for each
    dimension),

  - accounting for measurement error on the between- as well as
    within-level using multiple indicators (via specification of
    two-level confirmatory factor models),

  - all dynamic model parameters specified as constant or
    person-specific (random),

  - for bivariate VAR-models estimation of person-specific (random) log
    innovation covariances using a latent factor approach,

  - within-level latent interaction effects,

  - zero-inflation in time series variables via a censoring approach, as
    well as

  - between-level variables as predictors of random (individual)
    effects, and/or between-level outcomes predicted by random model
    parameters and additional between-level covariates.

- Generating data for all of the above modeling options (with additional
  convenience settings to study the models in simulation studies or for
  running power analyses based on specific parameter values).

- Accounting for non-equidistant measurements and overnight lags (by
  creating a time grid similar to the TINTERVAL option provided in
  Mplus).

- Running common posterior and model diagnostics analyses steps (as for
  objects of class stanfit) as well as specific posterior analyses for
  the mltsfit class.

- Multiple plotting options of the

  - specified model object as a path model,
  - group-level effects,
  - random effect parameters.

<div id="refs" class="references csl-bib-body hanging-indent"
data-entry-spacing="0" data-line-spacing="2">

<div id="ref-Asparouhov2018" class="csl-entry">

Asparouhov, T., Hamaker, E. L., & Muthén, B. (2018). Dynamic Structural
Equation Models. *Structural Equation Modeling: A Multidisciplinary
Journal*, *25*(3), 359–388.
<https://doi.org/10.1080/10705511.2017.1406803>

</div>

<div id="ref-rStan2023" class="csl-entry">

Stan Development Team. (2023a). *<span class="nocase">RStan: the R
interface to Stan</span>*. <https://mc-stan.org/>

</div>

<div id="ref-Stan2023" class="csl-entry">

Stan Development Team. (2023b). *<span class="nocase">Stan Modeling
Language Users Guide and Reference Manual, Version 2.31</span>*.
<https://mc-stan.org/>

</div>

</div>
