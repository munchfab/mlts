#' Fit Bayesian Multilevel Manifest or Latent Time-Series Models
#'
#' @param model data.frame. Output of \code{\link[mlts]{mlts_model}} and related functions.
#' @param data An object of class data.frame (or one that can be coerced to that
#' class) containing data of all variables used in the model.
#' @param id Character. The variable in `data` that identifies the person or observational
#' unit.
#' @param ts Character. The variable(s) in `data` that
#' contain the time-series construct. If multiple variables are provided in a
#' character vector, a vector autoregressive model is fit.
#' @param covariates Character. The covariate(s) in `data` used for prediction of
#' random effects.
#' @param outcomes Character. The outcome(s) in `data` that should be
#' predicted by random effects.
#' @param outcome_pred_btw Character.
#' @param center_covs Logical. Should covariates be centered before inclusion
#' in the model? Defaults to `TRUE`.
#' @param time Character. The variable in `data` that contains the (continuous) time.
#' @param tinterval The step interval for approximation for a continuous time
#' dynamic model. The smaller the step interval, the better the approximation.
#' @param beep tba.
#' @param days tba.
#' @param n_overnight_NAs tba.
#' @param na.rm tba.
#' @param iter A positive integer specifying the number of iterations for each
#' chain (including 50% used as warmup). The default is 500.
#' @param chains A positive integer specifying the number of Markov chains.
#' The default is 2.
#' @param cores The number of cores to use when executing the Markov chains in parallel.
#' The default is 2 (see \code{\link[rstan]{stan}}).
#' @param monitor_person_pars Logical. Should person parameters (i.e., values of the
#' latent variables) be stored? Default is FALSE.
#' @param print_message Logical. Print messages based on defined inputs (default = TRUE).
#' @param print_warning Logical. Print warnings based on defined inputs (default = TRUE).
#' @param fit_model Logical. Set to FALSE to avoid fitting the model which may be
#' helpful to inspect prepared data used for model estimation (default = TRUE).
#' @param ... Additional arguments passed to \code{\link[rstan]{sampling}}.
#'
#' @return An object of class `data.frame`.
#' @export
#'
#' @examples
#' \donttest{
#'  # build simple mlts model for two time-series variables
#'  model <- mlts_model(q = 2)
#'
#'  # fit model with (artificial) dataset ts_data
#'  fit <- mlts_fit(
#'    model = model,
#'    data = ts_data,
#'    ts = c("Y1", "Y2"), # time-series variables
#'    id = "ID", # identifier variable
#'    tinterval = 1 # interval for approximation of continuous-time dynamic model,
#'  )
#'
#'  # inspect model summary
#'  summary(fit)
#' }
#'
mlts_fit <- function(model,
                     data =NULL,
                     id,
                     ts,
                     covariates = NULL,
                     outcomes = NULL,
                     outcome_pred_btw = NULL,
                     center_covs = T,
                     time = NULL,
                     tinterval,
                     beep = NULL,
                     days = NULL,
                     n_overnight_NAs,
                     na.rm = F,
                     iter = 500,
                     chains = 2,
                     cores = 2,
                     monitor_person_pars = F,
                     fit_model = T,
                     print_message = T,
                     print_warning = T,
                     ...
){

  # eval the model
  infos <- mlts_model_eval(model)
  # Get the parameter table
  par_labels <- mlts_param_labels(model)

  # print information on indicators per dimension
  # print information on variables used for model estimation
  if(infos$isLatent == F){
    call_inds = c(
      "Time series variables as indicated by parameter subscripts: \n",
      unlist(lapply(1:infos$q, function(x){
        paste0("  ", x, " --> ", ts[x], "\n")
      }))
    )
  }
  if(infos$isLatent == T){
    call_inds = c(
      "Time series variables as indicated by parameter subscripts: \n",
      unlist(lapply(1:infos$q, function(x){
        paste0("  ", x, " --> ",
               paste0(ts[infos$indicators$q == x], collapse = " + "), "\n")
      }))
    )
  }
  cat(call_inds)

  # simulated data used
  data.simulated = ifelse(class(data)[1] == "mlts_simdata", TRUE, FALSE)

  # check if data is class "VARsimData"
  if(data.simulated == T){
    message("Simulated data provided:",
    "\nTrue scores used in the simulation will be added to the returned object.")

    par_labels <- merge(x = par_labels, data$model[,c("Param", "true.val")],
                       by = "Param", sort = F)

    # store true values of indivdual parameters
    re.trues <- data$RE.pars

    data <- data$data
    id <- "ID"
    beep <- "time"

    #### for now only use na.rm = T -version for simulated data
    na.rm <- TRUE
    }

  # Some initial checks:
  # avoiding specification of "covariates"- and "outcomes"-arguments,
  # if variable names in the model-object match the variables names in data
  if(is.null(covariates) & infos$n_cov>1){
    # check
    if(sum(!(infos$n_cov_vars %in% colnames(data)))>0){
      # not all variables found in the data
      stop("Not all between-level variables used as predictors of random effects
      in the model can be found in the data. You may need to rename the variable(s) in
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
      stop("Not all between-level outcome variables specified in the model
      can be found in the data. You may need to rename the variable(s) in
      the data or provide the name(s) via the `outcomes` argument.")

    } else {
      # create the necessary input to the outcomes-argument
      outcomes <- infos$out_var
      names(outcomes) <- infos$out_var
    }
  }

  if(is.null(outcome_pred_btw) & infos$n_z>0){
    # check
    if(sum(!(infos$n_z_vars %in% colnames(data)))>0){
      # not all variables found in the data
      stop("Not all additional between-level variables used in the outcome prediction
      model as specified in the model can be found in the data. You may need to rename the variable(s) in
      the data or provide the name(s) via the `outcome_pred_btw` argument.")

    } else {
      # create the necessary input to the outcomes-argument
      outcome_pred_btw <- infos$n_z_vars
      names(outcome_pred_btw) <- infos$n_z_vars
    }
  }

  # PREPARE DATA =========================================================
  data = prepare_data(data, id = id, ts = ts, time = time,
                      tinterval = tinterval, beep = beep,
                      days = days, n_overnight_NAs = n_overnight_NAs,
                      na.rm = na.rm, covariates = covariates,
                      outcomes = outcomes,
                      outcome_pred_btw = outcome_pred_btw)

  # ======================================================================

  # depending on model evaluate model type
  ## check if a measurement model is included
  isLatent <- ifelse(sum(model$Model == "Measurement")>0,TRUE,FALSE)


  # BY MODEL TYPE ========================================================
  # VAR(1) Models -------------------------------------------------------------

  ## Single-indicator VAR(1) model
  if(isLatent == F){
    # data preprocessing
    standata <- VARprepare(model = model, data = data, ts = ts,
                          covariates = covariates, outcomes = outcomes,
                          outcome_pred_btw = outcome_pred_btw,
                          center_covs = center_covs)

    # model fit
    pars <- c("gammas","b_fix", "sigma", "sd_R", "bcorr",
              "b_re_pred", "b_out_pred", "alpha_out", "sigma_out")
    if(monitor_person_pars == T){
      pars = c(pars, "b_free")
    }

    if(standata$n_inno_covs == 0 & standata$n_inno_cors == 0){
      if(fit_model==T){
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
    } else if(standata$n_inno_cors > 0){
      pars = c(pars,"bcorr_inn")
      if(fit_model==T){
        stanfit <- rstan::sampling(
          stanmodels$VAR_manifestCovsFix,
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
    } else if(standata$n_inno_covs > 0){
      if(fit_model==T){
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

  ## Multiple-indicator Model ==================================================
  if(isLatent == T){
    # data preprocessing
    standata <- VARprepare(model = model, data = data, ts = ts,
                           covariates = covariates, outcomes = outcomes,
                           outcome_pred_btw = outcome_pred_btw,
                           center_covs = center_covs)

    # model fit
    pars <- c("gammas","b_fix", "sigma", "sd_R", "bcorr",
              "b_re_pred", "b_out_pred", "alpha_out", "sigma_out",
              "alpha", "loadB", "sigmaB", "loadW", "sigmaW")
    if(monitor_person_pars == T){
      pars = c(pars, "b_free")
    }

    if(standata$n_inno_covs == 0 & standata$n_inno_cors == 0){
      if(fit_model==T){
        stanfit <- rstan::sampling(
          stanmodels$VAR_latent,
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
    } else if(standata$n_inno_cors > 0){
        pars = c(pars,"bcorr_inn")
        if(fit_model==T){
          stanfit <- rstan::sampling(
            stanmodels$VAR_latentCovsFix,
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
    } else if(standata$n_inno_covs > 0){
      if(fit_model==T){
        stanfit <- rstan::sampling(
          stanmodels$VAR_latentCovsRand,
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


  # ============================================================================


  if(fit_model == T){
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

    # add model parameter labels
    sums$Param_stan = row.names(sums)
    pop.sums <- merge(par_labels, y = sums, by = "Param_stan", sort = F)

    # create individual parameter summary table
    if(monitor_person_pars == TRUE){
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
        "num_id" = ID_new,
        "Param" = pars,
        sums.i
      )

      # if simulated data were used, add true values of person parameters
      if(data.simulated == T){
        sums.i.true = data.frame(
          "num_id" = rep(1:nrow(re.trues), ncol(re.trues)),
          "Param" = rep(colnames(re.trues), each = nrow(re.trues)),
          "true.val" = as.vector(re.trues)
        )
        sums.i = merge(x = sums.i.true, sums.i, by = c("num_id", "Param"), all = T)
      }


      # add original subject identifier
      id.match = unique(data[c(id, "num_id")])
      sums.i = merge(id.match, y = sums.i, by = "num_id", all = T)



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
  result <- list(
    "model" = model,
    "data" = data,
    "standata" = standata,
    "stanfit" = stanfit,
    "param.labels" = par_labels,
    "pop.pars.summary" = pop.sums,
    "posteriors" = posteriors,
    "person.pars.summary" = sums.i
  )

  # assign class
  class(result) <- "mltsfit"

  return(result)

}
