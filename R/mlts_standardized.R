#' Get Standardized Estimates for an mlts Model
#'
#' @param object `mltsfit`. Output of \code{\link[mlts]{mlts_model}} and related functions.
#' @param prob A value between 0 and 1 to indicate the width of the credible
#' interval. Default is .95.
#' @param digits Number of digits.
#' @param add_cluster_std logical. Include within-level standardized effects for each cluster
#' (defaults to `FALSE`).
#' @return An object of class `list`.
#' @export
#'
#' @examples
#' \donttest{
#' # build simple vector-autoregressive mlts model for two time-series variables
#' var_model <- mlts_model(q = 2)
#'
#' # fit model with (artificial) dataset ts_data
#' fit <- mlts_fit(
#'   model = var_model,
#'   data = ts_data,
#'   ts = c("Y1", "Y2"), # time-series variables
#'   id = "ID", # identifier variable
#'   tinterval = 1 # interval for approximation of continuous-time dynamic model,
#' )
#'
#' # inspect standardized parameter estimates
#' mlts_standardized(fit)
#' }
#'
mlts_standardized <- function(object, digit = 3, prob = .95, add_cluster_std = FALSE
){

  # make sure object is of class mltsfit
  if(class(object) != "mltsfit"){
    stop("Input of `object` should be of class 'mltsfit'.")
  }

  #
  if(length(object$person.pars.summary) == 1){
    stop("To obtain standardized estimates, set argument `monitor_person_pars = TRUE` in `mlts_fit()`")
  }

  # get model infos
  infos <- mlts_model_eval(object$model)

  # get CIs as specified by user
  alpha = 1 - prob
  probs = c(alpha/2, 1-alpha/2)
  prob.cols = paste0(c(100-(100-(alpha/2)*100), 100-(alpha/2)*100), "%")
  # result columns
  result.cols = c("Std.Est", "SD", prob.cols)

  # prepare results object
  results.std = list()
  btw.std = data.frame()
  within_std = data.frame()
  within_std_cluster = data.frame()

  # Between-level standardized effects =========================================
  ## Get posteriors of all random effects by parameter
  N = object$standata$N
  RE_par_SD = list()
  SDs = data.frame()
  chains <- as.numeric(object$stanfit@sim$chains)
  warmup <- as.numeric(object$stanfit@sim$warmup)
  iter <- as.numeric(object$stanfit@sim$iter) -warmup


  # # Old Version:
  # for(i in 1:object$standata$n_random){
  #   # create stan par names
  #   par_names = paste0("b_free[",1:N,",",i,"]")
  #   # get mcmc chains results for respective parameter
  #   poster = rstan::As.mcmc.list(object$stanfit, pars = par_names)
  #   # get individual parameter SD in each chain
  #   for(k in 1:chains){
  #     for(l in 1:nrow(poster[[k]])){
  #        SDs[l,i] = sd(poster[[k]][l,])
  #     }
  #     RE_par_SD[[k]] = SDs
  #   }
  # }

  # RE_par_SD_vec = c()
  # for(i in 1:object$standata$n_random){
  #   sds_vec = c()
  #   for(j in 1:chains){
  #     sds_vec = c(sds_vec, RE_par_SD[[j]][,i])
  #   }
  #   RE_par_SD_vec[i] = mean(sds_vec)
  # }
  #
  # RE_par_SD_vec


  ### NEW VERSION ::::::::

  RE_par_SD = c()
  for(i in 1:object$standata$n_random){
    # create stan par names
    par_names = paste0("b_free[",1:N,",",i,"]")
    poster = c()
    for(j in 1:chains){
      poster[j] = sd(rstan::get_posterior_mean(object = object$stanfit, pars = par_names)[,j])
    }
    RE_par_SD[i] = mean(poster)
  }

  ## random effect prediction
  if(object$standata$n_cov > 1){
    # prepare final object to store results
    re_pred_std = infos$RE.PREDS[, c("Type", "Param")]
    re_pred_std[,result.cols] = NA

    # start loop over all predictors
    for(i in 1:nrow(infos$RE.PREDS)){
     # get SD of covariate
     sd_x = sd(object$standata$W[,infos$RE.PREDS$pred_no[i]+1])
     # get position of RE on RE_par_SD
     sd_y_pos = infos$RE.PREDS$re_no[i]
     # get parameter label in stan model
     lab_x = object$pop.pars.summary$Param_stan[which(object$pop.pars.summary$Param == infos$RE.PREDS$Param[i])]
     # get unstandardized regression weight
     b_unstd = rstan::As.mcmc.list(object$stanfit, pars = lab_x)
     b_std_list = list()
     b_std = data.frame()

     for(k in 1:chains){
      # calculate standardized effect for all iterations chain-wise
      # formular: beta = b * (SD_x / SD_y)
      for(l in 1:nrow(b_unstd[[k]])){
        b_std[l,k] = b_unstd[[k]][l] * (sd_x / RE_par_SD[sd_y_pos])
      }
     }
     # calculate by iteration means
   #  b_std = apply(b_std, MARGIN = 1, FUN = mean)
     re_pred_std[i, result.cols] = round(c(
       mean(unlist(b_std)),
       sd(unlist(b_std)),
       quantile(unlist(b_std), c(probs))),digits = digit)
    }
    btw.std = rbind(btw.std, re_pred_std)
  }

  ## Outcome prediction
  if(object$standata$n_out > 0){
    # prepare final object to store results
    out_pred_std = infos$OUT[, c("Type", "Param")]
    out_pred_std[,result.cols] = NA

    # start loop over all predictors
    for(i in 1:nrow(infos$OUT)){
      # get SD of outcome
      sd_y = sd(object$standata$out[infos$OUT$out_var_no[i],])
      # get position of RE on RE_par_SD
      if(is.na(infos$OUT$Pred_Z[i])){
        sd_x = c()
        sd_x = RE_par_SD[infos$OUT$Pred_no[i]]
      } else {
        for(k in 1:chains){
          sd_x = sd(object$standata$Z[,infos$OUT$Pred_no[i]-infos$n_random])
        }
      }

      # get parameter label in stan model
      lab_x = object$pop.pars.summary$Param_stan[which(object$pop.pars.summary$Param == infos$OUT$Param[i])]
      # get unstandardized regression weight
      b_unstd = rstan::As.mcmc.list(object$stanfit, pars = lab_x)
      b_std_list = list()
      b_std = data.frame()

      for(k in 1:chains){
        # calculate standardized effect for all iterations chain-wise
        # formular: beta = b * (SD_x / SD_y)
        for(l in 1:nrow(b_unstd[[k]])){
          b_std[l,k] = b_unstd[[k]][l] * (sd_x / sd_y)
        }
      }
      # calculate by iteration means
      #b_std = apply(b_std, MARGIN = 1, FUN = mean)
      out_pred_std[i, result.cols] = round(c(
        mean(unlist(b_std)),
        sd(unlist(b_std)),
        quantile(unlist(b_std), c(probs))),digits = digit)
    }
    btw.std = rbind(btw.std, out_pred_std)
  }

  if(nrow(btw.std)>0){
    row.names(btw.std) <- NULL
    results.std[["Between-level standardized estimates"]] <- btw.std
  }













  # Average Within-Level Standardized Estimates of Dynamics ====================
  SD_y_id = data.frame()
  b_std = data.frame()
  within_std = infos$fix_pars_dyn[, c("Type", "Param")]
  std_individual <- array(dim = c(iter, chains, N))
  cluster_std <- list()
  for(p in 1:N){
    cluster_std[[p]] = infos$fix_pars_dyn[ ,c("Type", "Param")]
    }

  if(infos$isLatent == F){
    # first get individual SDs of time series variables
    for(i in 1:infos$q){
      for(pp in 1:N){
        SD_y_id[pp,i] = sd(object$data[object$data$num_id == pp, object$standata$ts[i]], na.rm = T)
      }
    }
    # calculate std estimates per person, averaged over chain and iteration
    for(j in 1:nrow(infos$fix_pars_dyn)){

      # get individual effect parameters
      if(infos$fix_pars_dyn$isRandom[j] == 1){
         re_par_no = infos$re_pars$par_no[infos$re_pars$Param == infos$fix_pars_dyn$Param[j]]
         param_stan = paste0("b_free[",1:N,",",re_par_no,"]")
      } else {
         fix_par_no = cumsum(infos$fix_pars_dyn$isRandom == 0)
         param_stan = paste0("b_fix[",fix_par_no[j],"]")
        }
      sd_y = SD_y_id[1:N,as.integer(infos$fix_pars_dyn$Dout[j])]
      sd_x = SD_y_id[1:N,as.integer(infos$fix_pars_dyn$Dpred[j])]
      b_unstd.list = rstan::As.mcmc.list(object$stanfit, pars = param_stan)
      b_std = data.frame()
      for(k in 1:chains){
        for(l in 1:iter){
          std_individual[l,k,1:N] = b_unstd.list[[k]][l,] * (sd_x /sd_y)
          b_std[l,k] = mean(std_individual[l,k,])
          }
        }
        #b_std = apply(b_std, MARGIN = 1, FUN = mean)
        within_std[j, result.cols] = round(c(
          mean(unlist(b_std)),
          sd(unlist(b_std)),
          quantile(unlist(b_std), c(probs))),digits = digit)
        for(p in 1:N){
          cluster_std[[p]][j,result.cols] = round(c(
            mean(unlist(std_individual[,,p])),
            sd(unlist(std_individual[,,p])),
            quantile(unlist(std_individual[,,p]), c(probs))),digits = digit)
        }
      }

  } else {  # for latent models ...

    print("Add standardized effects for latent model! ")
  }

  if(nrow(within_std)>0){
    row.names(within_std) <- NULL
    results.std[["Within-level standardidzed effects averaged over clusters"]] <- within_std
  }

  if(add_cluster_std==TRUE){
    results.std[["Standardized effects by cluster"]] <- cluster_std
  }

  return(results.std)

}
