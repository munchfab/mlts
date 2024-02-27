# mlts

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the development version of mlts from [GitHub](https://github.com/) with:

```r
# install.packages("devtools")
devtools::install_github("munchfab/mlts")
```

## To Do

### Model functions 
* `VARmodelBuild()`
  * [ ] rename to `mlts_model()`
  * [ ] include tests for unit testing with `testthat()`
  * [ ] write Vignettes
  * [ ] think about adding a model class argument to prepare future implementations of additional model types (e.g., TAR, modVAR, or rDSEM models)
* `VARmodelMeasurement()` 
  * [ ] rename to `mlts_model_measurement()`
  * [ ] include tests for unit testing with `testthat()`
  * [ ] write Vignettes
* `VARmodelBetween()`
  * [ ] rename to `mlts_model_betw()` 
  * include tests for unit testing with `testthat()`
* `VARmodelConstraints()` 
  * [ ] rename to `mlts_model_constraint()`
  * [ ] include tests for unit testing with `testthat()`
  * [ ] write Vignettes
* `VARmodelPriors()` 
  * [ ] rename to `mlts_model_priors()` 
  * [ ] include tests for unit testing with `testthat()`
  * [ ] write or include in Vignettes
* `VARmodelChecks()` 
  * [ ] rename to `mlts_model_check()`
  * [ ] function currently not integrated 
  * [ ] include more testing conditions where users could fail
  * [ ] include tests for unit testing with `testthat()`
* `VARmodelformula()` 
  * [ ] rename to `mlts_model_formula()`
  * [ ] include tests for unit testing with `testthat()`
  * [ ] include in Vignettes
* `VARmodelPaths()`
  * [ ] rename to `mlts_model_paths()` 
  * [ ] include tests for unit testing with `testthat()`
  * [ ] include in Vignettes

### Other Functions 
* `ARprepare()` and `VARprepare()`
  * `ARprepare()` is now obsolete --> all hendeled within `VARprepare()`
  * [ ] finish data preprocessing for AR- and VAR-models --> pass priors for constant model parameters to stan models (bfix)
  * [x] include `create_missings()` in data preparation function to allow for approximation of continuous time models
  * [ ] include tests for unit testing with `testthat()`
* `mlts_plot()`
  * [ ] think of a way to correctly print greek letters in plots   
  * [ ] include tests for unit testing with `testthat()`
* `VARfit()`
  * [ ] rename to `mlts_fit()` 
  * [ ] include tests for unit testing with `testthat()`
* `VARmodelParLabels()`
  * [ ] rename to `mlts_param_labels()`
  * [ ] should this be an explicit function (with help page) or just a background helper function? I prefer the latter
* `VARmodelSim()`:
  * [ ] rename to `mlts_sim()`  
  * [ ] include tests for unit testing with `testthat()`
  * [ ] write Vignettes
* `summary.mlts()`:
  * [ ] write a summary function to print a nice summary to console
  * possible options to include:
   * [ ] allow inclusion of standardized parameter estimates
 
### Stan Models 
* [ ] include user-specified priors on constant model parameters in all current models 
* [ ] bivarate VAR model with multiple-indicators and random innovation covaraince 
