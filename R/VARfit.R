#' Title
#'
#' @param VARmodel data.frame. Output of VARmodel-Functions.
#' @param data An object of class data.frame (or one that can be coerced to that
#' class) containing data of all variables used in the model.
#' @param ts.ind data.frame. Data input.
#' @param covariates data.frame. Data input.
#' @param outcomes data.frame. Data input.
#' @param outcome.pred.btw data.frame. Data input.
#' @param center.covs data.frame. Data input.
#' @param std.outcome data.frame. Data input.
#' @param iter A positive integer specifying the number of iterations for each
#' chain (including 50% used as warmup). The default is 500.
#' @param chains A positive integer specifying the number of Markov chains. The default is 2.
#' @param cores The number of cores to use when executing the Markov chains in parallel. The default is 2 (see \code{\link[rstan]{stan}}).
#' @param monitor.person.pars data.frame. Data input.
#' @param std.outcome data.frame. Data input.
#' @param printMessage logical. Print messages based on defined inputs (default = TRUE).
#' @param printWarning logical. Print messages based on defined inputs (default = TRUE).
#' @param fit.model logical. Set to FALSE to avoid fitting the model which may be
#' helpful to inspect prepared data used for model estimation (default = T).
#' @param ... Additional arguments passed to \code{\link[rstan]{sampling}}.
#'
#' @return An object of class `data.frame`.
#' @export
#'
VARfit <- function(VARmodel,
                   data =NULL,
                   ts.ind,
                   covariates = NULL,
                   outcomes = NULL,
                   outcome.pred.btw = NULL,
                   center.covs = T,
                   std.outcome = T,
                   iter = 500,
                   chains = 2,
                   cores = 2,
                   monitor.person.pars = F,
                   fit.model = T,
                   printMessage = T,
                   printWarning = T,
                   ...
){

  # eval the model
  infos <- VARmodelEval(VARmodel)
  # Get the parameter table
  par_labels <- VARmodelParLabels(VARmodel)


  # check if data is class "VARsimData"
  if(class(data) == "VARsimData"){
    message("Simulated data provided: True scores used in the simulation will
            added to the returned object.")

    par_labels <- merge(x = par_labels, data$VARmodel[,c("Param", "true.val")],
                       by = "Param", sort = F)
    data = data$data
    }

  # Some initial checks:
  # avoiding specification of "covariates"- and "outcomes"-arguments,
  # if variable names in the VARmodel-object match the variables names in data
  if(is.null(covariates) & infos$n_cov>1){
    # check
    if(sum(!(infos$n_cov_vars %in% colnames(data)))>0){
      # not all variables found in the data
      stop("Not all between-level variables used as predictors of random effects
      in the VARmodel can be found in the data. You may need to rename the variable(s) in
      the data or provide the name(s) via the `covariates` argument.")

    } else {
      # create the necessary input to the outcomes-argument
      covariates <- infos$n_cov_vars
      names(covariates) <- infos$n_cov_vars
    }
  }

  if(is.null(outcomes) & infos$n_out>0){
    # check
    if(sum(!(infos$out_var %in% colnames(data)))>0){
      # not all variables found in the data
      stop("Not all between-level outcome variables specified in the VARmodel
      can be found in the data. You may need to rename the variable(s) in
      the data or provide the name(s) via the `outcomes` argument.")

    } else {
      # create the necessary input to the outcomes-argument
      outcomes <- infos$out_var
      names(outcomes) <- infos$out_var
    }
  }

  if(is.null(outcome.pred.btw) & infos$n_z>0){
    # check
    if(sum(!(infos$n_z_vars %in% colnames(data)))>0){
      # not all variables found in the data
      stop("Not all additional between-level variables used in the outcome prediction
      model as specified in the VARmodel can be found in the data. You may need to rename the variable(s) in
      the data or provide the name(s) via the `outcome.pred.btw` argument.")

    } else {
      # create the necessary input to the outcomes-argument
      outcome.pred.btw <- infos$n_z_vars
      names(outcome.pred.btw) <- infos$n_z_vars
    }
  }


  # depending on VARmodel evaluate model type
  ## number of time-varying constructs
  isAR <- sum(VARmodel$Param_Label == "Trait" & VARmodel$Type == "Fix effect") == 1

  ## number if random effects
  n_random <- sum(VARmodel$isRandom, na.rm = T)

  ## check if a measurement model is included
  isLatent <- ifelse(sum(VARmodel$Model == "Measurement")>0,TRUE,FALSE)


  # BY MODEL TYPE ========================================================

  # AR(1) Models -------------------------------------------------------------
  ## Single-indicator AR(1) model with at least two random effects
  if(isAR == T & n_random > 1 & isLatent == F){
    # data preprocessing
    if(printMessage==T){
      standata <- ARprepare(VARmodel = VARmodel, data = data, ts.ind = ts.ind,
                           covariates = covariates, outcomes = outcomes,
                           outcome.pred.btw = outcome.pred.btw,
                           center.covs = center.covs, std.outcome = std.outcome)
    } else {
      standata <- suppressMessages(
        ARprepare(VARmodel = VARmodel, data = data, ts.ind = ts.ind,
                  covariates = covariates, outcomes = outcomes,
                  outcome.pred.btw = outcome.pred.btw,
                  center.covs = center.covs, std.outcome = std.outcome))
    }

    # parameter to monitor
    pars <- c("gammas", "sigma", "sd_R", "bcorr", "b_re_pred", "b_out_pred", "alpha_out", "sigma_out")
    if(monitor.person.pars == T){
      pars <- c(pars, "b_free")
    }

    # fit model
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
      standata <- ARprepare(VARmodel = VARmodel, data = data, ts.ind = ts.ind,
                           covariates = covariates, outcomes = outcomes,
                           outcome.pred.btw = outcome.pred.btw,
                           center.covs = center.covs, std.outcome = std.outcome)
    } else {
      standata <- suppressMessages(
        ARprepare(VARmodel = VARmodel, data = data, ts.ind = ts.ind,
                  covariates = covariates, outcomes = outcomes,
                  outcome.pred.btw = outcome.pred.btw,
                  center.covs = center.covs, std.outcome = std.outcome))
    }

    # model fit
    pars <- c("gammas", "ar", "sigma", "sd_R", "b_re_pred", "b_out_pred", "alpha_out", "sigma_out")
    if(monitor.person.pars == T){
      pars <- c(pars, "b_free")
    }

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
    standata <- data

    # model fit
  }

  ## Multiple-indicator AR(1) model with random intercepts only
  if(isAR == T & n_random == 1 & isLatent == T){
    # data preprocessing
    standata <- data

    # model fit
  }

  # VAR(1) Models -------------------------------------------------------------
  ## Single-indicator VAR(1) model
  if(isAR == F & isLatent == F){
    # data preprocessing
    standata <- VARprepare(VARmodel = VARmodel, data = data, ts.ind = ts.ind,
                          covariates = covariates, outcomes = outcomes,
#                          outcome.pred.btw = outcome.pred.btw,
                          center.covs = center.covs, std.outcome = std.outcome)

    # model fit
    pars <- c("gammas", "b_fix", "btw_free", "sd_R","sigma", "bcorr")
    if(monitor.person.pars == T){
      pars = c(pars, "b_free")
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
        stanfit <- NULL
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
        stanfit <- NULL
      }
    }


  }

  ## Multiple-indicator VAR(1) model
  if(isAR == F & isLatent == T){
    # data preprocessing
    standata <- data

    # model fit
  }


  if(fit.model == T){
    # add posteriors with adapted names
    posteriors <- rstan::extract(stanfit, inc_warmup = F, permuted = F,
                             pars = par_labels$Param_stan)
    dimnames(posteriors)$parameters <- par_labels$Param

    # create a summary table using the monitor-function in rstan
    sums <- rstan::monitor(stanfit, print = F)

    # get a subset of outputs
    cols <- c("mean", "sd", "2.5%", "50%", "97.5%", "n_eff",
             "Rhat", "Bulk_ESS", "Tail_ESS")
    sums <- round(sums[1:dim(sums)[1], cols],3)

    # add VARmodel parameter labels
    sums$Param_stan = row.names(sums)
    pop.sums <- merge(par_labels, y = sums, by = "Param_stan", sort = F)

    # create individual parameter summary table
    if(monitor.person.pars == TRUE){
      sums.i = sums[startsWith(sums$Param_stan, "b_free"),1:ncol(sums)]
      # extract infos
      pars <- gsub(sums.i$Param_stan, pattern = "b_free[", replacement = "", fixed = T)
      pars <- gsub(pars, pattern = "]", replacement = "", fixed = T)
      ID_new <- sapply(pars, function(x){strsplit(x,split = ",")[[1]][1]})
      pars <- sapply(pars, function(x){as.integer(strsplit(x,split = ",")[[1]][2])})
      pars <- sapply(pars, function(x){infos$re_pars$Param[x]})

      row.names(sums.i) <- NULL
      sums.i$Param_stan <- NULL
      sums.i = cbind(
        "ID_new" = ID_new,
        "Param" = pars,
        sums.i
      )
    } else {
        # we could think about adding the posterior means here
        sums.i <- NA
    }

  } else {
    posteriors <- NA
    pop.sums <- NA
    sums.i <- NA
  }


  # combine preprocessed data and fitted stan object in a list object
  VARresult <- list(
    "VARmodel" = VARmodel,
    "standata" = standata,
    "stanfit" = stanfit,
    "param.labels" = par_labels,
    "pop.pars.summary" = pop.sums,
    "posteriors" = posteriors,
    "person.pars.summary" = sums.i
  )

  return(VARresult)

}
