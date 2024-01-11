#' Title
#'
#' @param VARmodel data.frame. Output of VARmodel-Functions.
#' @param data data.frame. Data input.
#' @param printMessage logical. Print messages based on defined inputs (default = TRUE).
#' @param fit.model logical. Set to FALSE to avoid fitting the model which may be
#' helpful to inspect prepared data used for model estimation (default = T).
#' @param ... Additional arguments passed to rstan `sampling`-function.
#'
#' @return An object of class `data.frame`.
#' @export
VARfit <- function(VARmodel, data, printMessage = T, printWarning = T,
                   ts.ind, covariates = NULL, outcomes = NULL,
                   center.covs = T, std.outcome = T, iter = 500, chains = 2, cores = 2,
                   monitor.person.pars = F,
                   fit.model = T,
                   ...
){

  # depending on VARmodel evaluate model type
  ## number of time-varying constructs
  isAR = sum(VARmodel$Param_Label == "Trait" & VARmodel$Type == "Fix effect") == 1

  ## number if random effects
  n_random = sum(VARmodel$isRandom, na.rm = T)

  ## check if a measurement model is included
  isLatent = ifelse(sum(VARmodel$Model == "Measurement")>0,TRUE,FALSE)


  # BY MODEL TYPE ========================================================

  # AR(1) Models -------------------------------------------------------------
  ## Single-indicator AR(1) model with at least two random effects
  if(isAR == T & n_random > 1 & isLatent == F){
    # data preprocessing
    if(printMessage==T){
      standata = ARprepare(VARmodel = VARmodel, data = data, ts.ind = ts.ind,
                           covariates = covariates, outcomes = outcomes,
                           center.covs = center.covs, std.outcome = std.outcome)
    } else {
      standata = suppressMessages(
        ARprepare(VARmodel = VARmodel, data = data, ts.ind = ts.ind,
                  covariates = covariates, outcomes = outcomes,
                  center.covs = center.covs, std.outcome = std.outcome))
    }

    # model fit
    pars <- c("btw_pred", "sigma", "bcorr")

    stanfit <- rstan::sampling(
      stanmodels$AR_manifest,
      data = standata,
      pars = pars,
      iter = iter,
      cores = cores,
      chains = chains,
      ...
    )
  }

  ## Single-indicator AR(1) model with random intercepts only
  if(isAR == T & n_random == 1 & isLatent == F){
    # data preprocessing
    if(printMessage==T){
      standata = ARprepare(VARmodel = VARmodel, data = data, ts.ind = ts.ind,
                           covariates = covariates, outcomes = outcomes,
                           center.covs = center.covs, std.outcome = std.outcome)
    } else {
      standata = suppressMessages(
        ARprepare(VARmodel = VARmodel, data = data, ts.ind = ts.ind,
                  covariates = covariates, outcomes = outcomes,
                  center.covs = center.covs, std.outcome = std.outcome))
    }

    # model fit
    pars <- c("btw_pred", "ar", "sd_noise", "sd_R")

    stanfit <- rstan::sampling(
      stanmodels$AR_manifest_intOnly,
      data = standata,
      pars = pars,
      iter = iter,
      cores = cores,
      chains = chains,
      ...
    )
  }

  ## Multiple-indicator AR(1) model with at least two random effects
  if(isAR == T & n_random > 1 & isLatent == T){
    # data preprocessing
    standata = data

    # model fit
  }

  # ## Multiple-indicator AR(1) model with random intercepts only
  # if(isAR == T & n_random == 1 & isLatent == T){
  #   # data preprocessing
  #   standata = data
  #
  #   # model fit
  }

  # VAR(1) Models -------------------------------------------------------------
  ## Single-indicator VAR(1) model
  if(isAR == F & isLatent == F){
    # data preprocessing
    standata = VARprepare(VARmodel = VARmodel, data = data, ts.ind = ts.ind,
                          covariates = covariates, outcomes = outcomes,
                          center.covs = center.covs, std.outcome = std.outcome)

    # model fit
    pars <- c("gammas", "b_fix","btw_free", "sd_R","sigma", "bcorr")
    if(monitor.person.pars == T){
      pars <- c(pars, "btw_free")
    }

    if(standata$n_inno_covs == 0){
      if(fit.model==T){
        stanfit <- rstan::sampling(
          stanmodels$VAR_manifest,
          data = standata,
          pars = pars,
          iter = iter,
          cores = cores,
          chains = chains,
          ...
        )
      } else {
        stanfit = NULL
      }
    } else {
      if(fit.model==T){
        stanfit <- rstan::sampling(
          stanmodels$VAR_manifestCovsRand,
          data = standata,
          pars = pars,
          iter = iter,
          cores = cores,
          chains = chains,
          ...
        )
      } else {
        stanfit = NULL
      }
    }


  }

  ## Multiple-indicator VAR(1) model
  if(isAR == F & isLatent == T){
    # data preprocessing
    standata = data

    # model fit
  }


  # combine preprocessed data and fitted stan object in a list object
  VARresult = list(
    "VARmodel" = VARmodel,
    "standata" = standata,
    "stanfit" = stanfit
  )

  return(VARresult)

}
