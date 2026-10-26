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
mlts_standardize_within <- function(object, digits = 3, prob = .95, add_cluster_std = FALSE,
                                    has_cluster_samples = TRUE
){

  # get model infos
  infos <- mlts_model_eval(object$model)

  # get CIs as specified by user
  alpha = 1 - prob
  probs = c(alpha/2, 1-alpha/2)
  prob.cols = paste0(c(100-(100-(alpha/2)*100), 100-(alpha/2)*100), "%")
  # result columns
  if(infos$G == 1){
    result.cols = c("Std.Est", "SD", prob.cols)
  } else {
    result.cols = c("Group", "Std.Est", "SD", prob.cols)
  }


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

  if(has_cluster_samples == FALSE){
    chains = 1
    iter = 1
  }

  # Average Within-Level Standardized Estimates of Dynamics ====================
  n_dim = infos$q
  within_std = list()
  intermediate_cols = c("isRandom","isAR","isINT","isRDSEM","Dout","Dpred","no")
  for(gg in 1:infos$G){
    within_std[[gg]] = infos$fix_pars_dyn[, c("Type", "Param",intermediate_cols)]

    if(nrow(infos$fix_pars_rdsem)>0){
      within_std[[gg]] <- rbind(
        infos$fix_pars_rdsem[, c("Type", "Param",intermediate_cols)],
        infos$fix_pars_dyn  [, c("Type", "Param",intermediate_cols)])
    }
  }

  cluster_std <- list()
  for(p in 1:N){
    if(nrow(infos$fix_pars_rdsem)>0){
      cluster_std[[p]] = rbind(
        infos$fix_pars_rdsem[, c("Type", "Param",intermediate_cols)],
        infos$fix_pars_dyn[, c("Type", "Param",intermediate_cols)])
    } else {
      cluster_std[[p]] = infos$fix_pars_dyn[, c("Type", "Param",intermediate_cols)]
    }
  }

  ### get variances ============================================================
  SD_y_id = std_get_within_vars(infos = infos, object = object, has_cluster_samples = has_cluster_samples)

  # ADD LATENT VERSION!
  if(infos$n_rdsem.pars > 0){
    SD_yresid_id = std_get_within_vars_resid(infos = infos, object = object, has_cluster_samples = has_cluster_samples)
  }

  ### calculate standardized effects ===========================================

  # create a list with posterior draws of std effect for each cluster
  bstd_all = lapply(1:nrow(within_std[[1]]), function(j){
    b = array(dim = c(N, iter*chains))
    b_std = array(dim = c(N, iter*chains))
    if(has_cluster_samples == FALSE){
      means = rstan::get_posterior_mean(object$stanfit)
    }

    for(p in 1:N){

      # pick parameter
      if(within_std[[1]]$isRandom[j] == 1){
        par_no = within_std[[1]]$no[j]
        param_stan = paste0("b_free[",p,",",par_no,"]")
      } else {
        par_no = cumsum(within_std[[1]]$isRandom == 0)[j]
        g_id = object$standata$g_id[p]
        param_stan = paste0("b_fix[",g_id,",",par_no[j],"]")
      }
      # extract unstandardized estimates
      if(has_cluster_samples == TRUE){
        b[p,] = rstan::extract(object$stanfit, pars = param_stan)[[1]]
      } else {
        b[p,] = mean(means[rownames(means) == param_stan,])
      }

      # pick SDs
      dv_num = as.integer(within_std[[1]]$Dout[j])
      iv_num = as.integer(within_std[[1]]$Dpred[j])
      SDy = SD_y_id[dv_num,p,]
      SDx = SD_y_id[iv_num,p,]

      # check adjust in case of residual DSEM path
      if(dv_num %in% infos$fix_pars_rdsem$Dout & within_std[[1]]$isRDSEM[j] != 1){
        SDy = SD_yresid_id[dv_num,p,]
      }
      if(iv_num %in% infos$fix_pars_rdsem$Dout & within_std[[1]]$isRDSEM[j] != 1){
        SDx = SD_yresid_id[iv_num,p,]
      }

      # calculate standardized estimate
      b_std[p,] = b[p,] * SDx / SDy

      # interaction effects
      if(within_std[[1]]$isINT[j] == 1){
        iv2_num = infos$fix_pars_dyn$Dpred2[infos$fix_pars_dyn$Param == within_std[[1]]$Param[j]]
        SDx2 = SD_y_id[iv2_num,p,]

        if(iv2_num %in% infos$fix_pars_rdsem$Dout & within_std[[1]]$isRDSEM[j] != 1){
          SDx = SD_yresid_id[iv2_num,p,]
        }

        # calculate standardized estimate
        b_std[p,] = b[p,] / SDy * SDx * SDx2
      }
    }

    return(b_std)

  })

  names(bstd_all) <- within_std[[1]]$Param



  # get posterior draws of average within-cluster standardized effects
  # (by group if necessary)
  g_id = object$standata$g_id
  for(gg in 1:infos$G){

    within_std[[gg]] = do.call(rbind, lapply(1:length(bstd_all), function(x){
      b_std_group = bstd_all[[x]][object$standata$g_id == gg,]
      if(has_cluster_samples == TRUE){
        b_std_average = apply(b_std_group, MARGIN = 2, FUN = mean)
      } else {
        b_std_average = mean(b_std_group)
      }

      out = cbind.data.frame(
        "Param" = names(bstd_all)[x],
        "Group" = object$standata$group_lab[gg],
        round(mean(b_std_average), digits = digits),
        round(stats::sd(b_std_average), digits = digits),
        t(round(stats::quantile(b_std_average, c(probs)),digits = digits)))

      if(has_cluster_samples == FALSE){
        out[,4:6] <- NA
      }

      if(infos$G == 1){
        out$Group <- NULL
      }
      colnames(out)[2:ncol(out)] <- result.cols

      return(out)
    }))

  }

#  get cluster-specific estimates and summary statistics
for(p in 1:N){
  for(j in 1:nrow(cluster_std[[p]])){
    g_id = object$standata$g_id[p]
    if(infos$G > 1){
      cluster_std[[p]][j,"Group"] = g_id
    }
    cluster_std[[p]][j,result.cols[result.cols!="Group"]] = round(c(
              mean(bstd_all[[j]][p,]),
              stats::sd(bstd_all[[j]][p,]),
              stats::quantile(bstd_all[[j]][p,], c(probs))),digits = digits
              )
  }
  cluster_std[[p]][,intermediate_cols] <- NULL
}

  # bring list into data frame format
  if(infos$G == 1){
    within_std <- within_std[[1]]
  } else {
    within_std <- do.call(rbind, within_std)
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


# mlts_standardize_within <- function(object, digits = 3, prob = .95, add_cluster_std = FALSE
# ){
#
#   # get model infos
#   infos <- mlts_model_eval(object$model)
#
#   # get CIs as specified by user
#   alpha = 1 - prob
#   probs = c(alpha/2, 1-alpha/2)
#   prob.cols = paste0(c(100-(100-(alpha/2)*100), 100-(alpha/2)*100), "%")
#   # result columns
#   if(infos$G == 1){
#     result.cols = c("Std.Est", "SD", prob.cols)
#   } else {
#     result.cols = c("Group", "Std.Est", "SD", prob.cols)
#   }
#
#
#   # prepare results object
#   results.std = list()
#   btw.std = data.frame()
#   within_std = data.frame()
#   within_std_cluster = data.frame()
#
#   # get information on fitted stan model
#   N = object$standata$N
#   chains <- as.numeric(object$stanfit@sim$chains)
#   warmup <- as.numeric(object$stanfit@sim$warmup)
#   iter <- as.numeric(object$stanfit@sim$iter)-warmup
#
#
#   # Average Within-Level Standardized Estimates of Dynamics ====================
#   n_dim = infos$q
#   SD_y_id = array(dim = c(n_dim, N, (chains*iter)))
#   within_std = list()
#   for(gg in 1:infos$G){
#     within_std[[gg]] <- infos$fix_pars_dyn[, c("Type", "Param")]
#   }
#
#   cluster_std <- list()
#   for(p in 1:N){
#     cluster_std[[p]] = infos$fix_pars_dyn[ ,c("Type", "Param")]
#   }
#
#   if(infos$isLatent == FALSE){
#     # first get individual SDs of time series variables
#     for(i in 1:infos$q){
#       for(pp in 1:N){
#         SD_y_id[i,pp,] = stats::sd(object$data[object$data$num_id == pp, object$standata$ts[i]], na.rm = TRUE)
#       }
#     }
#     # calculate std estimates per person
#     for(j in 1:nrow(infos$fix_pars_dyn)){
#       b = array(dim = c(N, iter*chains))
#       b_std = array(dim = c(N, iter*chains))
#       # get individual effect parameters
#       if(infos$fix_pars_dyn$isRandom[j] == 1){
#         re_par_no = infos$re_pars$par_no[infos$re_pars$Param == infos$fix_pars_dyn$Param[j]]
#         for(p in 1:N){
#           param_stan = paste0("b_free[",p,",",re_par_no,"]")
#           b[p,] = rstan::extract(object$stanfit, pars = param_stan)[[1]]
#           b_std[p,] = b[p,] *
#             SD_y_id[as.integer(infos$fix_pars_dyn$Dpred[j]),p,] / # sd_x
#             SD_y_id[as.integer(infos$fix_pars_dyn$Dout[j]),p,]   # sd_y
#           if(infos$fix_pars_dyn$isINT[j] == 1){  # interaction effects
#             b_std[p,] = b[p,] /
#               SD_y_id[as.integer(infos$fix_pars_dyn$Dout[j]),p,] *  # sd_y
#               SD_y_id[as.integer(infos$fix_pars_dyn$Dpred[j]),p,] *  # sd_x1
#               SD_y_id[as.integer(infos$fix_pars_dyn$Dpred2[j]),p,]   # sd_x2
#             }
#           }
#         } else {
#           fix_par_no = cumsum(infos$fix_pars_dyn$isRandom == 0)
#           for(p in 1:N){
#             # get fixed effect - depending on group id
#             g_id = object$standata$g_id[p]
#             param_stan = paste0("b_fix[",g_id,",",fix_par_no[j],"]")
#
#             b[p,] = rstan::extract(object$stanfit, pars = param_stan)[[1]]
#             b_std[p,] = b[p,] *
#               SD_y_id[as.integer(infos$fix_pars_dyn$Dpred[j]),p,] / # sd_x
#               SD_y_id[as.integer(infos$fix_pars_dyn$Dout[j]),p,]   # sd_y
#
#             if(infos$fix_pars_dyn$isINT[j] == 1){  # interaction effects
#               b_std[p,] = b[p,] /
#                 SD_y_id[as.integer(infos$fix_pars_dyn$Dout[j]),p,] *  # sd_y
#                 SD_y_id[as.integer(infos$fix_pars_dyn$Dpred[j]),p,] *  # sd_x1
#                 SD_y_id[as.integer(infos$fix_pars_dyn$Dpred2[j]),p,]   # sd_x2
#               }
#           }
#         }
#       # calculate average standardized effect per iteration
#       if(infos$G == 1){
#         b_std_average = apply(b_std, MARGIN = 2, FUN = mean)
#         within_std[[1]][j, result.cols] = round(c(
#           mean(b_std_average),
#           stats::sd(b_std_average),
#           stats::quantile(b_std_average, c(probs))),digits = digits)
#       } else {
#         for(gg in 1:infos$G){
#           b_std_group = b_std[object$standata$g_id == gg,]
#           b_std_average = apply(b_std_group, MARGIN = 2, FUN = mean)
#           within_std[[gg]][j, result.cols] = round(c(
#             object$standata$group_lab[gg],
#             mean(b_std_average),
#             stats::sd(b_std_average),
#             stats::quantile(b_std_average, c(probs))),digits = digits)
#           }
#         }
#
#
#
#       # get cluster-specific estimates
#       if(infos$G == 1){
#         for(p in 1:N){
#           cluster_std[[p]][j,result.cols] = round(c(
#             mean(b_std[p,]),
#             stats::sd(b_std[p,]),
#             stats::quantile(b_std[p,], c(probs))),digits = digits)
#             }
#       } else {
#         for(p in 1:N){
#           g_id = object$standata$g_id[p]
#           cluster_std[[p]][j,result.cols] = round(c(
#             object$standata$group_lab[g_id],
#             mean(b_std[p,]),
#             stats::sd(b_std[p,]),
#             stats::quantile(b_std[p,], c(probs))),digits = digits)
#           }
#         }
#
#       }
#
#   } else if(object$standata$standardized == 0){  # check if SDs of latent variables are available
#     warning("Variance of latent factor scores not available for standardization
#             of dynamic model parameters. Refit the model with get_SD_latent = TRUE
#             to obtain standardized estimates.")
#
#     } else {  # run standardization
#     k_etaW_index = 1 # index latent factor constructs
#     SD_y_id = array(dim=c(n_dim, N, (chains*iter)))
#
#     # first get SDs of time series variables
#     for(i in 1:infos$q){
#       if(object$standata$D_np[i] == 1){  # for constructs with single-indicator
#         q_pos_ts = object$standata$D_pos_is_SI[i]
#         for(pp in 1:N){
#           SD_y_id[i,pp,] = stats::sd(object$data[object$data$num_id == pp, object$standata$ts[q_pos_ts]], na.rm = TRUE)
#         }
#       } else {  # constructs with multiple indicators
#         for(pp in 1:N){
#           # get sds of latent variables
#           etaW_sd_lab = paste0("SD_etaW_i[",k_etaW_index,",",pp,"]")
#           SD_y_id[i,pp,] = rstan::extract(object$stanfit, pars = etaW_sd_lab)[[1]]
#         }
#         # update index
#         k_etaW_index = k_etaW_index + 1
#       }
#     }
#
#     # calculate std estimates per person, averaged over chain and iteration
#     # calculate std estimates per person
#     for(j in 1:nrow(infos$fix_pars_dyn)){
#       b = array(dim = c(N, iter*chains))
#       b_std = array(dim = c(N, iter*chains))
#       # get individual effect parameters
#       if(infos$fix_pars_dyn$isRandom[j] == 1){
#         re_par_no = infos$re_pars$par_no[infos$re_pars$Param == infos$fix_pars_dyn$Param[j]]
#         for(p in 1:N){
#           param_stan = paste0("b_free[",p,",",re_par_no,"]")
#           b[p,] = rstan::extract(object$stanfit, pars = param_stan)[[1]]
#           b_std[p,] = b[p,] *
#             SD_y_id[as.integer(infos$fix_pars_dyn$Dpred[j]),p,] / # sd_x
#             SD_y_id[as.integer(infos$fix_pars_dyn$Dout[j]),p,]   # sd_y
#           if(infos$fix_pars_dyn$isINT[j] == 1){  # interaction effects
#             b_std[p,] = b[p,] /
#               SD_y_id[as.integer(infos$fix_pars_dyn$Dout[j]),p,] *  # sd_y
#               SD_y_id[as.integer(infos$fix_pars_dyn$Dpred[j]),p,] *  # sd_x1
#               SD_y_id[as.integer(infos$fix_pars_dyn$Dpred2[j]),p,]   # sd_x2
#           }
#         }
#       } else {
#         fix_par_no = cumsum(infos$fix_pars_dyn$isRandom == 0)
#         param_stan = paste0("b_fix[",fix_par_no[j],"]")
#         b[,] = rstan::extract(object$stanfit, pars = param_stan)[[1]]
#         b_std[p,] = b[p,] *
#           SD_y_id[as.integer(infos$fix_pars_dyn$Dpred[j]),p,] / # sd_x
#           SD_y_id[as.integer(infos$fix_pars_dyn$Dout[j]),p,]   # sd_y
#         if(infos$fix_pars_dyn$isINT[j] == 1){  # interaction effects
#           b_std[p,] = b[p,] /
#             SD_y_id[as.integer(infos$fix_pars_dyn$Dout[j]),p,] *  # sd_y
#             SD_y_id[as.integer(infos$fix_pars_dyn$Dpred[j]),p,] *  # sd_x1
#             SD_y_id[as.integer(infos$fix_pars_dyn$Dpred2[j]),p,]   # sd_x2
#         }
#       }
#       # calculate average standardized effect per iteration
#       b_std_average = apply(b_std, MARGIN = 2, FUN = mean)
#       within_std[[1]][j, result.cols] = round(c(
#         mean(b_std_average),
#         stats::sd(b_std_average),
#         stats::quantile(b_std_average, c(probs))),digits = digits)
#
#       # get cluster-specific estimates
#       for(p in 1:N){
#         cluster_std[[p]][j,result.cols] = round(c(
#           mean(b_std[p,]),
#           stats::sd(b_std[p,]),
#           stats::quantile(b_std[p,], c(probs))),digits = digits)
#       }
#     }
#   }
#
#   # bring list into data frame format
#   if(infos$G == 1){
#     within_std <- within_std[[1]]
#   } else {
#     within_std <- do.call(rbind, within_std)
#   }
#
#   if(nrow(within_std)>0){
#     row.names(within_std) <- NULL
#     results.std[["Within-level standardidzed effects averaged over clusters"]] <- within_std
#   }
#
#   if(add_cluster_std==TRUE){
#     results.std[["Within-level standardized effects by cluster"]] <- cluster_std
#   }
#
#   return(results.std)
#
# }
