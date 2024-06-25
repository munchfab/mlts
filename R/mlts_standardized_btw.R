#' Helper function for standardizing within-level parameters
#'
#' @param object `mltsfit`. Output of \code{\link[mlts]{mlts_model}} and related functions.
#' @param digits Number of digits. Default is 3.
#' @param prob A value between 0 and 1 to indicate the width of the credible
#' interval. Default is .95.
#'
#' @return An object of class `list` containing between-level standardized estimates.
#' @noRd
#'
mlts_standardized_btw <- function(object, digits = 3, prob = .95
){

  # make sure object is of class mltsfit
  # if(class(object) != "mltsfit")
  if (!inherits(object, "mltsfit")) {
    stop("Input of `object` should be of class 'mltsfit'.")
  }

  # get model infos
  infos <- mlts_model_eval(object$model)

  # get CIs as specified by user
  alpha = 1 - prob
  probs = c(alpha/2, 1-alpha/2)
  prob.cols = paste0(c(100-(100-(alpha/2)*100), 100-(alpha/2)*100), "%")
  # result columns
  result.cols = c("Std.Est", "SD", prob.cols)
  result.list <- list() # final object
  result <- data.frame()

  # Between-level standardized effects =========================================
  ## Get posteriors of all random effects by parameter
  chains <- as.numeric(object$stanfit@sim$chains)
  warmup <- as.numeric(object$stanfit@sim$warmup)
  iter <- as.numeric(object$stanfit@sim$iter) -warmup
  n_random <- object$standata$n_random
  re_sds = array(dim = c(chains, iter, n_random))
  N = object$standata$N
  # get SDs of RE predictor variables
  Wvars_sds = apply(object$standata$W, MARGIN = 2, FUN = stats::sd)

  # get/calculate model-implied RE variances ----
  ## get labels of (residual) random effect SDs
  labs_re = paste0("sd_R[",1:n_random,"]")
  sigma_RE = rstan::extract(object$stanfit, pars = labs_re)
  # number of covariates used for prediction of random effects
  n_cov = object$standata$n_cov
  # as default use sigma of random effects and update respectively
  Var_RE = sigma_RE
  for(i in 1:n_random){
    Var_RE[[i]] = Var_RE[[i]]^2
  }

  # Part 1: Standardized estimates for random effects regressed on covariate(s)
  # Get variance of Y predicted based on unstandardized regression weights
  if(object$standata$n_cov > 1){
    # loop over random effects
    for(i in 1:n_random){
      # get number and position of covariates
      n_cov_pred = sum(object$standata$n_cov_mat[,2] == i)
      n_cov_pos = object$standata$n_cov_mat[which(object$standata$n_cov_mat[,2] == i),1]

      if(n_cov_pred > 0){
        # store y predicted
        y_pred = array(data = 0, dim = c(chains*iter, N))
        for(j in 1:n_cov_pred){
          cov_mat = object$standata$n_cov_mat
          b_pos = which(cov_mat[,1] == n_cov_pos[j] & cov_mat[,2] == i)
          b <- rstan::extract(object$stanfit, pars = paste0("b_re_pred[",b_pos,"]"))
          for(k in 1:(iter*chains)){
            y_pred[k,] = y_pred[k,] + object$standata$W[,n_cov_pos[j]] * b[[1]][k]
          }
        }
        # get variance of predicted scores in each iteration
        y_pred_var = array(dim = chains*iter)
        for(k in 1:(iter*chains)){
          y_pred_var[k] = stats::var(y_pred[k,])
        }

        # add variance of predicted scores to sigma in each iteration
        Var_RE[[i]] = Var_RE[[i]] + y_pred_var
      }
    }


    # now use SDs of RE to standardize regression parameters:
    # prepare object to store results
    re_pred_std = infos$RE.PREDS[, c("Type", "Param")]
    re_pred_std[,result.cols] = NA
    for(k in 1:nrow(infos$RE.PREDS)){

      # get SD of covariate
      sd_x = Wvars_sds[object$standata$n_cov_mat[k,1]]
      # get position of RE on RE_par_SD
      sd_y = sqrt(Var_RE[[object$standata$n_cov_mat[k,2]]])
      # get parameter label in stan model
      b <- rstan::extract(object$stanfit, pars = paste0("b_re_pred[",k,"]"))
      b_std <- b[[1]] * sd_x / sd_y
      # save summary statistics
      re_pred_std[k, result.cols] = round(c(
        mean(unlist(b_std)),
        stats::sd(unlist(b_std)),
        stats::quantile(unlist(b_std), c(probs))),digits = digits)
    }

    result = rbind(result, re_pred_std)
  }

  # Part 2: Standardized estimates for outcomes regressed on random effects
  if(object$standata$n_out > 0){
    # prepare final object to store results
    out_pred_std = infos$OUT[, c("Type", "Param")]
    out_pred_std[,result.cols] = NA

    # get SDs of outcomes
    SDs_out = apply(object$standata$out, FUN = stats::sd, MARGIN = 1)
    # add variances of additional covariates used as predictor
    if(object$standata$n_z > 0){
      for(i in 1:object$standata$n_z){
        Var_RE[[n_random+i]] <- stats::var(object$standata$Z[,i])
      }
    }

    for(i in 1:nrow(infos$OUT)){
      # get unstandardized estimate
      lab_out = object$pop.pars.summary$Param_stan[which(object$pop.pars.summary$Param == infos$OUT$Param[i])]
      b = rstan::extract(object$stanfit, pars = c(lab_out))
      sd_y = SDs_out[infos$OUT$out_var_no[i]]
      sd_x = sqrt(Var_RE[[infos$OUT$Pred_no[i]]])
      b_std = b[[1]] * sd_x / sd_y
      # save summary statistics
      out_pred_std[i, result.cols] = round(c(
        mean(unlist(b_std)),
        stats::sd(unlist(b_std)),
        stats::quantile(unlist(b_std), c(probs))),digits = digits)
    }

    row.names(result) <- NULL
    result = rbind(result, out_pred_std)
  }

  ## Standardization of constant dynamic parameters ----------------------------
  # check if single- or multiple-indicator model
  isLatent <- ifelse(sum(object$model$Model == "Measurement")>0,TRUE,FALSE)

  # run standardization for single-indicator models ----------------------------
  if(isLatent == FALSE){
  # check that for each dimension:
  # all dynamic parameters are fixed
  # all innovation variances of dependent and independent dimensions are fixed
  fix_dyn = infos$fix_pars_dyn[infos$fix_pars_dyn$isRandom == 0,]
  if(nrow(fix_dyn) > 0){
    # use average of observed intraindividual variance for standardization of
    # constant dynamic parameters
    VarYw = array(dim = c(infos$q))
    for(i in 1:infos$q){
      ivar = c()
      for(p in 1:N){
        ivar[p] = stats::var(object$data[object$data$num_id==p,object$standata$ts[i]], na.rm=TRUE)
      }
      VarYw[i] = mean(ivar)
    }


    # run checks looped over dynamic parameters
    for(i in 1:nrow(fix_dyn)){
      # 1. check if all dynamic parameters predicting the respective construct
      # are constant across subjects
      dim_out = fix_dyn$Dout[i]
      all_dynPars_fixed = sum(infos$fix_pars_dyn[infos$fix_pars_dyn$Dout == dim_out,"isRandom"]) == 0
      # 2. check if innovation variances of involved constructs are constant
      dim_pred = infos$fix_pars_dyn$Dpred[infos$fix_pars_dyn$Dout == dim_out]
      # check based on parameter labels
      InnoVar_labs = paste0("sigma_",dim_pred)
      all_InnoVars_fixed = sum(!(InnoVar_labs %in% infos$fix_pars$Param)) == 0

      # run standardization
      if(all_dynPars_fixed & all_InnoVars_fixed){
        # get unstandardized estimates per iteration
        par_label = fix_dyn$Param[i]
        par_stan = object$param.labels$Param_stan[object$param.labels$Param == par_label]
        b = rstan::extract(object$stanfit, pars = c(par_stan))
        sd_x = sqrt(VarYw[as.integer(fix_dyn$Dpred[i])])
        sd_y = sqrt(VarYw[as.integer(fix_dyn$Dout[i])])
        b_std <- b[[1]] * sd_x / sd_y
        b_std = round(c(
          mean(unlist(b_std)),
          stats::sd(unlist(b_std)),
          stats::quantile(unlist(b_std), c(probs))),digits = digits)
        b_std = cbind.data.frame(
          "Type" = "Dynamic",
          "Param" = par_label,
          t(b_std)
        )
        colnames(b_std) = c("Type", "Param", result.cols)
        result = rbind(result, b_std)
      }
    }
  }
  }


  # run standardization for multiple-indicator models --------------------------
  if(isLatent == TRUE & object$standata$standardized == 0){
    warning("Variance of latent factor scores not available for standardization
            of dynamic model parameters. Refit the model with get_SD_latent = TRUE
            to obtain standardized estimates.")
  } else if(isLatent == TRUE & object$standata$standardized == 1){

    # check that for each dimension:
    # all dynamic parameters are fixed
    # all innovation variances of dependent and independent dimensions are fixed
    fix_dyn = infos$fix_pars_dyn[infos$fix_pars_dyn$isRandom == 0,]
    k_etaW_index = 1 # index latent factor constructs
    if(nrow(fix_dyn) > 0){
      VarYw = array(dim = c(infos$q, (chains*iter)))

      for(i in 1:infos$q){
        # if number of indiactors is 1, use average intraindividual variance
        if(object$standata$D_np[i] == 1){
          ivar = c()
          for(p in 1:N){
            ivar[p] = stats::var(object$data[object$data$num_id==p,object$standata$ts[i]], na.rm=TRUE)
          }
          VarYw[i,] = mean(ivar)
        } else {
          # get sds of latent variables
          etaW_sd = rstan::extract(object$stanfit, pars = paste0("SD_etaW[",k_etaW_index,"]"))
          VarYw[i,] = etaW_sd[[1]]^2
          # update index
          k_etaW_index = k_etaW_index + 1
        }
      }

      # run checks looped over dynamic parameters
      for(i in 1:nrow(fix_dyn)){
        # 1. check if all dynamic parameters predicting the respective construct
        # are constant across subjects
        dim_out = fix_dyn$Dout[i]
        all_dynPars_fixed = sum(infos$fix_pars_dyn[infos$fix_pars_dyn$Dout == dim_out,"isRandom"]) == 0
        # 2. check if innovation variances of involved constructs are constant
        dim_pred = infos$fix_pars_dyn$Dpred[infos$fix_pars_dyn$Dout == dim_out]
        # check based on parameter labels
        InnoVar_labs = paste0("sigma_",dim_pred)
        all_InnoVars_fixed = sum(!(InnoVar_labs %in% infos$fix_pars$Param)) == 0

        # run standardization
        if(all_dynPars_fixed & all_InnoVars_fixed){
          # get unstandardized estimates per iteration
          par_label = fix_dyn$Param[i]
          par_stan = object$param.labels$Param_stan[object$param.labels$Param == par_label]
          b = rstan::extract(object$stanfit, pars = c(par_stan))
          sd_x = sqrt(VarYw[as.integer(fix_dyn$Dpred[i])])
          sd_y = sqrt(VarYw[as.integer(fix_dyn$Dout[i])])
          b_std <- b[[1]] * sd_x / sd_y
          b_std = round(c(
            mean(unlist(b_std)),
            stats::sd(unlist(b_std)),
            stats::quantile(unlist(b_std), c(probs))),digits = digits)
          b_std = cbind.data.frame(
            "Type" = "Dynamic",
            "Param" = par_label,
            t(b_std)
          )
          colnames(b_std) = c("Type", "Param", result.cols)
          result = rbind(result, b_std)
        }
      }
    }
  }

  return(result)
}
