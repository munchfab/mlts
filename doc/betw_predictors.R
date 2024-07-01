## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(mlts)

## -----------------------------------------------------------------------------
load("gait_data.rda")
head(gait_data)

## -----------------------------------------------------------------------------
ar1_model1 <- mlts_model(q = 1)

## -----------------------------------------------------------------------------
ar1_model1

## ----eval=FALSE---------------------------------------------------------------
#  mlts_model_paths(ar1_model1)

## ----echo=FALSE, out.width="75%", fig.align='center'--------------------------
knitr::include_graphics(c("../vignettes/pathmodel_ar1.png"))

## ----eval=FALSE---------------------------------------------------------------
#  ar1_fit1 <- mlts_fit(
#    model = ar1_model1,
#    data = gait_data,
#    id = "subject",
#    ts = "stride_interval",
#    monitor_person_pars = TRUE,
#    iter = 1000
#  )

## ----echo=FALSE---------------------------------------------------------------
load("../vignettes/ar1_fit1.rda")

## ----eval=FALSE---------------------------------------------------------------
#  summary(ar1_fit1)

## ----echo=FALSE---------------------------------------------------------------
cat(readLines("../vignettes/betw_preds_fit1.txt"), sep = "\n")

## -----------------------------------------------------------------------------
ar1_model2 <- mlts_model(q = 1, ranef_pred = "healthy")

## -----------------------------------------------------------------------------
ar1_model2

## ----eval=FALSE---------------------------------------------------------------
#  mlts_model_paths(ar1_model2)

## ----echo=FALSE, out.width="30%", fig.align='center'--------------------------
knitr::include_graphics(c("../vignettes/pathmodel_betw_preds.png"))

## ----eval=FALSE---------------------------------------------------------------
#  ar1_fit2 <- mlts_fit(
#    model = ar1_model2,
#    data = gait_data,
#    id = "subject",
#    ts = "stride_interval",
#    center_covs = FALSE,
#    monitor_person_pars = TRUE,
#    iter = 1000
#  )

## ----echo=FALSE---------------------------------------------------------------
load("../vignettes/ar1_fit2.rda")

## ----eval=FALSE---------------------------------------------------------------
#  summary(ar1_fit2)

## ----echo=FALSE---------------------------------------------------------------
cat(readLines("../vignettes/betw_preds_fit2.txt"), sep = "\n")

## ----echo=FALSE---------------------------------------------------------------
load("../vignettes/ar1_fit2_std.rda")

## ----eval=FALSE---------------------------------------------------------------
#  mlts_standardized(ar1_fit2, what = "both")

## ----echo=FALSE---------------------------------------------------------------
ar1_fit2_std

