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

* `ARprepare()` and `VARprepare()` 
  [x] finish data preprocessing for AR- and VAR-models
  * include `create_missings()` in data preparation function to allow for approximation of continuous time models
  * include tests for unit testing with `testthat()`
* `mlts_plot()` 
  * include tests for unit testing with `testthat()`
* `VARfit()` 
  * include tests for unit testing with `testthat()`
* `VARmodelBetween()` 
  * include tests for unit testing with `testthat()`
* `VARmodelBuild()` 
  * include tests for unit testing with `testthat()`
  * write Vignettes
* `VARmodelChecks()` 
  * include more testing conditions where users could fail
  * include tests for unit testing with `testthat()`
* `VARmodelConstraints()` 
  * include tests for unit testing with `testthat()`
  * write Vignettes
* `VARmodelformula()` 
  * include tests for unit testing with `testthat()`
  * include in Vignettes
* `VARmodelMeasurement()` 
  * include tests for unit testing with `testthat()`
  * write Vignettes
* `VARmodelParLabels()`
* `VARmodelPaths()` 
  * include tests for unit testing with `testthat()`
  * include in Vignettes
* `VARmodelPriors()` 
  * include tests for unit testing with `testthat()`
  * write or include in Vignettes
* `VARmodelSim()`: 
  * include tests for unit testing with `testthat()`
  * write Vignettes
* Decide on package name and rename functions accordingly?
