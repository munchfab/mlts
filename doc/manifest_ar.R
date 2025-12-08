## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup--------------------------------------------------------------------
library(mlts)

## -----------------------------------------------------------------------------
ar1_model <- mlts_model(q = 1)

## -----------------------------------------------------------------------------
ar1_model

## ----eval=FALSE---------------------------------------------------------------
# mlts_model_formula(ar1_model)

## ----echo=FALSE, out.width="50%", fig.align='center', fig.width=7, fig.height=5----
mlts_paths(ar1_model)

## ----echo=FALSE---------------------------------------------------------------
load("../data/ar1_data.rda")

## -----------------------------------------------------------------------------
head(ar1_data)

## ----eval=FALSE---------------------------------------------------------------
# ar1_fit <- mlts_fit(
#   model = ar1_model,
#   data = ar1_data,
#   id = "ID",
#   ts = "Y1",
#   iter = 4000
# )
# 
# saveRDS(ar1_fit, file = "../vignettes/ar1_fit.rds")

## ----echo=FALSE---------------------------------------------------------------
ar1_fit <- readRDS("../vignettes/ar1_fit.rds")

## -----------------------------------------------------------------------------
summary(ar1_fit)

