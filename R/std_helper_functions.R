


std_get_within_vars <- function(infos, object, has_cluster_samples = TRUE){
  N = object$standata$N
  chains <- as.numeric(object$stanfit@sim$chains)
  warmup <- as.numeric(object$stanfit@sim$warmup)
  iter <- as.numeric(object$stanfit@sim$iter)-warmup
  n_dim = infos$q

  if(has_cluster_samples == FALSE){
    chains = 1
    iter = 1
  }

  SD_y_id = array(dim = c(n_dim, N, (chains*iter)))

  if(infos$isLatent == FALSE){
    # first get individual SDs of time series variables
    for(i in 1:infos$q){
      for(pp in 1:N){
        SD_y_id[i,pp,] = stats::sd(object$data[object$data$num_id == pp, object$standata$ts[i]], na.rm = TRUE)
      }
    }
  }

  if(infos$isLatent == TRUE & object$standata$standardized == 0){
    warning("Variance of latent factor scores not available for standardization
            of dynamic model parameters. Refit the model with get_SD_latent = TRUE
            to obtain standardized estimates.")
  }

  if(infos$isLatent == TRUE & object$standata$standardized == 1){
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
      ##### ACHTUNG:
      ##### Für rDSEM-Modelle müsste hier noch die durch Prädiktor aufgeklärte
      ##### Varianz hinzugefügt werden.
    }

  }

  return(SD_y_id)
  }


std_get_within_vars_resid<- function(infos, object, has_cluster_samples = TRUE){

  N = object$standata$N
  chains <- as.numeric(object$stanfit@sim$chains)
  warmup <- as.numeric(object$stanfit@sim$warmup)
  iter <- as.numeric(object$stanfit@sim$iter)-warmup

  if(has_cluster_samples == FALSE){
    iter = 1
    chains = 1
    means = rstan::get_posterior_mean(object$stanfit)
  }

  SD_y_id = std_get_within_vars(infos = infos, object = object, has_cluster_samples = has_cluster_samples)
  SD_yresid_id = SD_y_id

  if(infos$isLatent == FALSE){

    for(pp in 1:N){
      # extract observations
      obs_id = nrow(object$data[object$data$num_id == pp,])
      y_obs = array(dim = c(infos$q, obs_id))
      for(i in 1:infos$q){
        y_obs[i,] = object$data[object$data$num_id == pp, object$standata$ts[i]]
      }

      y_resid = lapply(1:iter, function(x){
        y_obs
        })
      rdsem_paths = infos$fix_pars_rdsem
      rdsem_paths$fix_par_no = cumsum(rdsem_paths$isRandom == 0)

      for(j in 1:nrow(rdsem_paths)){
        iv_num = as.integer(rdsem_paths$Dpred[j])
        dv_num = as.integer(rdsem_paths$Dout[j])

        # get regression parameter:
        if(rdsem_paths$isRandom[j] == 1){
          par_no = rdsem_paths$no[j]
          param_stan = paste0("b_free[",pp,",",par_no,"]")
        } else {
          g_id = object$standata$g_id[pp]
          par_no = rdsem_paths$fix_par_no[j]
          param_stan = paste0("b_fix[",g_id,",",par_no,"]")
        }

        if(has_cluster_samples == TRUE){
          b = rstan::extract(object$stanfit, pars = param_stan)[[1]]
        } else {
          b = mean(means[rownames(means) == param_stan,])
        }

        for(ii in 1:iter){
          iv     = y_resid[[ii]][iv_num,]
          dv     = y_resid[[ii]][dv_num,]
          y_resid[[ii]][dv_num,] =  y_resid[[ii]][dv_num,] - b[ii]*y_resid[[ii]][iv_num,]
        }
      }

      # calculate SDs of residual variables
      for(j in 1:infos$q){
        for(ii in 1:iter){
          SD_yresid_id[j,pp,ii] = stats::sd(y_resid[[ii]][j,], na.rm = T)
        }
      }
    }
  }




  # if(infos$isLatent == TRUE){
  #
  #   ###### TBA -> NEEDS TO BE EXTRACTED FROM STAN-MODEL!
  #
  # }


  return(SD_yresid_id)

}

