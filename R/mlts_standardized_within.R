#' Title
#'
#' @param object `mltsfit`. Output of \code{\link[mlts]{mlts_model}} and related functions.
#' @param digits Number of digits. Default is 3.
#' @param prob A value between 0 and 1 to indicate the width of the credible
#' interval. Default is .95.
#' @param add_cluster_std logical. If `what = "within"`, within-level standardized effects for each cluster
#' are included in the output (defaults to `FALSE`).
#'
#' @return An object of class `list` containing within-level standardized estimates.
#' @noRd
#'
mlts_standardize_within <- function(object, digits = 3, prob = .95, add_cluster_std = FALSE
){

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

  # get information on fitted stan model
  N = object$standata$N
  chains <- as.numeric(object$stanfit@sim$chains)
  warmup <- as.numeric(object$stanfit@sim$warmup)
  iter <- as.numeric(object$stanfit@sim$iter)-warmup


  # Average Within-Level Standardized Estimates of Dynamics ====================
  n_dim = infos$q
  SD_y_id = array(dim = c(n_dim, N, (chains*iter)))
  within_std = infos$fix_pars_dyn[, c("Type", "Param")]
  cluster_std <- list()
  for(p in 1:N){
    cluster_std[[p]] = infos$fix_pars_dyn[ ,c("Type", "Param")]
  }

  if(infos$isLatent == FALSE){
    # first get individual SDs of time series variables
    for(i in 1:infos$q){
      for(pp in 1:N){
        SD_y_id[i,pp,] = stats::sd(object$data[object$data$num_id == pp, object$standata$ts[i]], na.rm = TRUE)
      }
    }
    # calculate std estimates per person
    for(j in 1:nrow(infos$fix_pars_dyn)){
      b = array(dim = c(N, iter*chains))
      b_std = array(dim = c(N, iter*chains))
      # get individual effect parameters
      if(infos$fix_pars_dyn$isRandom[j] == 1){
        re_par_no = infos$re_pars$par_no[infos$re_pars$Param == infos$fix_pars_dyn$Param[j]]
        for(p in 1:N){
          param_stan = paste0("b_free[",p,",",re_par_no,"]")
          b[p,] = rstan::extract(object$stanfit, pars = param_stan)[[1]]
          b_std[p,] = b[p,] *
            SD_y_id[as.integer(infos$fix_pars_dyn$Dpred[j]),p,] / # sd_x
            SD_y_id[as.integer(infos$fix_pars_dyn$Dout[j]),p,]   # sd_y
          }
        } else {
          fix_par_no = cumsum(infos$fix_pars_dyn$isRandom == 0)
          param_stan = paste0("b_fix[",fix_par_no[j],"]")
          b[,] = rstan::extract(object$stanfit, pars = param_stan)[[1]]
          b_std[p,] = b[p,] *
            SD_y_id[as.integer(infos$fix_pars_dyn$Dpred[j]),p,] / # sd_x
            SD_y_id[as.integer(infos$fix_pars_dyn$Dout[j]),p,]   # sd_y
        }
      # calculate average standardized effect per iteration
      b_std_average = apply(b_std, MARGIN = 2, FUN = mean)
      within_std[j, result.cols] = round(c(
        mean(b_std_average),
        stats::sd(b_std_average),
        stats::quantile(b_std_average, c(probs))),digits = digits)

      # get cluster-specific estimates
      for(p in 1:N){
        cluster_std[[p]][j,result.cols] = round(c(
          mean(b_std[p,]),
          stats::sd(b_std[p,]),
          stats::quantile(b_std[p,], c(probs))),digits = digits)
      }
    }

  } else if(object$standata$standardized == 0){  # check if SDs of latent variables are available
    warning("Variance of latent factor scores not available for standardization
            of dynamic model parameters. Refit the model with get_SD_latent = TRUE
            to obtain standardized estimates.")

    } else {  # run standardization
    k_etaW_index = 1 # index latent factor constructs
    SD_y_id = array(dim=c(n_dim, N, (chains*iter)))

    # first get SDs of time series variables
    for(i in 1:infos$q){
      if(object$standata$D_np[i] == 1){  # for constructs with single-indicator
        q_pos_ts = object$standata$D_pos_is_SI[i]
        for(pp in 1:N){
          SD_y_id[i,pp,] = stats::sd(object$data[object$data$num_id == pp, object$standata$ts[q_pos_ts]], na.rm = TRUE)
        }
      } else {  # constructs with multiple indicators
        for(pp in 1:N){
          # get sds of latent variables
          etaW_sd_lab = paste0("SD_etaW_i[",k_etaW_index,",",pp,"]")
          SD_y_id[i,pp,] = rstan::extract(object$stanfit, pars = etaW_sd_lab)[[1]]
        }
        # update index
        k_etaW_index = k_etaW_index + 1
      }
    }

    # calculate std estimates per person, averaged over chain and iteration
    # calculate std estimates per person
    for(j in 1:nrow(infos$fix_pars_dyn)){
      b = array(dim = c(N, iter*chains))
      b_std = array(dim = c(N, iter*chains))
      # get individual effect parameters
      if(infos$fix_pars_dyn$isRandom[j] == 1){
        re_par_no = infos$re_pars$par_no[infos$re_pars$Param == infos$fix_pars_dyn$Param[j]]
        for(p in 1:N){
          param_stan = paste0("b_free[",p,",",re_par_no,"]")
          b[p,] = rstan::extract(object$stanfit, pars = param_stan)[[1]]
          b_std[p,] = b[p,] *
            SD_y_id[as.integer(infos$fix_pars_dyn$Dpred[j]),p,] / # sd_x
            SD_y_id[as.integer(infos$fix_pars_dyn$Dout[j]),p,]   # sd_y
        }
      } else {
        fix_par_no = cumsum(infos$fix_pars_dyn$isRandom == 0)
        param_stan = paste0("b_fix[",fix_par_no[j],"]")
        b[,] = rstan::extract(object$stanfit, pars = param_stan)[[1]]
        b_std[p,] = b[p,] *
          SD_y_id[as.integer(infos$fix_pars_dyn$Dpred[j]),p,] / # sd_x
          SD_y_id[as.integer(infos$fix_pars_dyn$Dout[j]),p,]   # sd_y
      }
      # calculate average standardized effect per iteration
      b_std_average = apply(b_std, MARGIN = 2, FUN = mean)
      within_std[j, result.cols] = round(c(
        mean(b_std_average),
        stats::sd(b_std_average),
        stats::quantile(b_std_average, c(probs))),digits = digits)

      # get cluster-specific estimates
      for(p in 1:N){
        cluster_std[[p]][j,result.cols] = round(c(
          mean(b_std[p,]),
          stats::sd(b_std[p,]),
          stats::quantile(b_std[p,], c(probs))),digits = digits)
      }
    }
  }

  if(nrow(within_std)>0){
    row.names(within_std) <- NULL
    results.std[["Within-level standardidzed effects averaged over clusters"]] <- within_std
  }

  if(add_cluster_std==TRUE){
    results.std[["Within-level standardized effects by cluster"]] <- cluster_std
  }

  return(results.std)

}
