mlts_standardize_within <- function(object, digit = 3, prob = .95, add_cluster_std = FALSE
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

  if(infos$isLatent == F){
    # first get individual SDs of time series variables
    for(i in 1:infos$q){
      for(pp in 1:N){
        SD_y_id[i,pp,] = sd(object$data[object$data$num_id == pp, object$standata$ts[i]], na.rm = T)
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
        sd(b_std_average),
        quantile(b_std_average, c(probs))),digits = digit)

      # get cluster-specific estimates
      for(p in 1:N){
        cluster_std[[p]][j,result.cols] = round(c(
          mean(b_std[p,]),
          sd(b_std[p,]),
          quantile(b_std[p,], c(probs))),digits = digit)
      }
    }

  } else if(object$standata$standardized == 0){  # check if SDs of latent variables are available
  } else {  # run standardization
    k_etaW_index = 1 # index latent factor constructs
    SD_y_id = array(dim=c())

    # first get SDs of time series variables
    for(i in 1:infos$q){
      for(pp in 1:N){
        SD_y_id[pp,i] = sd(object$data[object$data$num_id == pp, object$standata$ts[i]], na.rm = T)
      }
    }
    # calculate std estimates per person, averaged over chain and iteration
    for(j in 1:nrow(infos$fix_pars_dyn)){

      # check


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
