mlts_sim_within <- function(
    infos,
    burn.in,
    N,
    group_ids,
    TP,       # single integer or a vector of length N
    btw,
    mm_pars = NULL,
    cor_pars = NULL,
    exogenous = NULL){

  # turn TP into vector
  if(length(TP) == 1){TP <- rep(TP, times = N)}
  names(TP) = 1:N

  # prepare data frame
  within = data.frame(
    "ID"   = unlist(lapply(1:N, function(x){rep(x,TP[x])})),
    "time" = unlist(lapply(unname(TP), function(x){1:x}))
  )

  if(infos$G > 1){
    within$group = unlist(lapply(1:N, function(x){rep(group_ids[x],TP[x])}))
  }


  q = infos$q                # number of constructs
  y_cols = paste0("Y",1:q)   # prepare columns
  within[, y_cols] = NA

  # get positions to fill transition matrix
  dyn = infos$fix_pars_dyn
  if(nrow(dyn)>0){
    dyn$tran_pos = NA
    dyn$tran_pos = as.integer(dyn$Dpred) + q*(as.numeric(dyn$Lag)-1)
  }
  dyn_int = dyn[dyn$isINT == 1,]
  dyn = dyn[dyn$isINT != 1,]
  # remove t0-effects at this points
  dyn_t0 = dyn[dyn$Lag == 0,]
  dyn = dyn[dyn$Lag!=0,]

  # check w_cen vars
  wcen_logical = infos$is_wcen == 1
  n_wcen = sum(wcen_logical)
  exo_pos = 1

  for(i in 1:N){
    NT <- TP[i] + burn.in
    use_group <- group_ids[i]

    # build person-specific transition matrix
    transition = matrix(nrow = q, ncol = q*infos$maxLag, data = 0)
    y = matrix(nrow = NT, ncol = q)

    # build person-specific prediction error matrix
    innoVars.i = btw[i,infos$innos_pos]

    for(d in 1:q){ # loop over number of constructs
      if(infos$is_wcen[d] == 1){
        # parameter positions on btw-matrix
        par_pos = dyn$no[dyn$Dout==d]
        # fill transition matrix with person-specific parameter values
        transition[d,as.integer(dyn$tran_pos[dyn$Dout==d])] = btw[i,par_pos]

        if(infos$innos_rand[infos$D_cen_pos[d]] == 1){
          innoVars.i[infos$D_cen_pos[d]] = exp(innoVars.i[infos$D_cen_pos[d]]) # retransform log innovations
        } else{
          innoVars.i[infos$D_cen_pos[d]] = innoVars.i[infos$D_cen_pos[d]]^2
        }
      }
    }

    if(q == 1 | sum(wcen_logical)==1){
      inno_var_mat = matrix(data = innoVars.i, nrow = 1, ncol = 1)
    } else {
      inno_var_mat = diag(innoVars.i)

      if(infos$n_inno_cors > 0){
        for(xx in 1:q){
          for(yy in 1:q){
            par_is_there = cor_pars$sample[cor_pars$group==use_group & cor_pars$Param == paste0("r.zeta_",xx,yy)]
            if(xx < yy & length(par_is_there) > 0){
              x_pos = infos$D_cen_pos[xx]
              y_pos = infos$D_cen_pos[yy]
              cor = cor_pars$sample[cor_pars$group==use_group & cor_pars$Param == paste0("r.zeta_",xx,yy)]
              cov = cor * sqrt(inno_var_mat[x_pos,x_pos]) * sqrt(inno_var_mat[y_pos,y_pos])
              inno_var_mat[x_pos,y_pos] <- inno_var_mat[y_pos,x_pos] <- cov
            }
          }
        }
      }
    }

    # generate the within-level process in a loop over time points
    for(t in 1:NT){

      if(t <= infos$maxLag | t == 1){
        # starting values
        init = matrix(data = NA, nrow = 1, ncol = q)
        init[,wcen_logical] = mvtnorm::rmvnorm(n = 1, mean = rep(0, n_wcen), sigma = inno_var_mat)
        if(any(infos$is_wcen == 0)){
          init[,!wcen_logical] = exogenous[exo_pos, ]
        }
        y[t,] = init

        # start the process when sufficient initial values are generated
      } else {
        y_lag = c()
        y_lag = as.vector(y[t-1,1:q])            # create a lagged vector

        if(infos$maxLag>1){                      # ... extend for higher-order
          for(ll in 2:infos$maxLag){
            y_lag = c(y_lag, as.vector(y[t-ll,1:q]))
          }
        }

        if(q > 1 | infos$maxLag > 1){
          # get expected values using matrix multiplication
          if( ncol(transition)  == 0 ) {
            y[t,1:q] = 0
          } else {
            y[t,1:q] = y_lag %*% t(transition)
          }
          if(any(infos$is_wcen == 0)){
            y[t,!wcen_logical] = exogenous[exo_pos, ]
          }

        } else if( ncol(transition) >0 ) {
          # or for AR(1):
          y[t,] = y_lag * transition
        } else if( ncol(transition) == 0) {
          # for location scale models
          y[t,] = 0
        }

        # add innovations
        y[t,wcen_logical] = y[t,wcen_logical] + mvtnorm::rmvnorm(n = 1, mean = rep(0, n_wcen), sigma = inno_var_mat)

        # add contemporaneous effects
        if(nrow(dyn_t0) != 0){
          for(k in 1:nrow(dyn_t0)){
            dv <- as.integer(dyn_t0$Dout[k])
            iv <- as.integer(dyn_t0$Dpred[k])
            y[t,dv] = y[t,dv] + btw[i,dyn_t0$no[k]] * y[t,iv]
          }
        }

        # add interaction effects
        if(nrow(dyn_int) > 0){
          for(k in 1:nrow(dyn_int)){
            dv <- as.integer(dyn_int$Dout[k])
            iv1 <- as.integer(dyn_int$Dpred[k])
            iv2 <- as.integer(dyn_int$Dpred2[k])
            lag1 <- as.integer(dyn_int$Lag[k])
            lag2 <- as.integer(dyn_int$Lag2[k])
            y[t,dv] = y[t,dv] + btw[i,dyn_int$no[k]] * y[t-lag1,iv1] * y[t-lag2,iv2]
          }
        }

        # for bivariate VAR-models with random innovation covariance factor:
        if(infos$q >= 2 & infos$n_inno_covs == 1){
          inno_t = stats::rnorm(n = 1, mean = 0, sd = sqrt(exp(btw[i,infos$inno_cov_pos])))
          y[t,1:2] = y[t,1:2] + infos$inno_cov_load * inno_t
        }
      }
      exo_pos = exo_pos +1
    }

    # add structural regressions paths (rDSEM)
    if(infos$n_rdsem.pars > 0){
      for(l in 1:nrow(infos$fix_pars_rdsem)){
        dv <- as.integer(infos$fix_pars_rdsem$Dout[l])
        iv <- as.integer(infos$fix_pars_rdsem$Dpred[l])
        y[,dv] = y[,dv] + btw[i,infos$fix_pars_rdsem$no[l]] * y[,iv]
      }
    }

    # remove burn-in
    within[within$ID==i, y_cols] = y[(burn.in+1) : (burn.in+TP[i]),]


    # add trait scores (for manifest indicators)
    if(infos$isLatent == FALSE){
      for(j in 1:infos$q){
        if(infos$is_wcen[j] == TRUE){
          within[within$ID==i,y_cols[j]] = btw[i,infos$D_cen_pos[j]] + within[within$ID==i,y_cols[j]]
        }
      }
    }
  }

  # --------


  ##### add measurement model here -------------------------------------------
  # create manifest indicator scores
  if(infos$isLatent == TRUE){
    N_inds = max(infos$indicators$p_pos)
    for(i in 1:N_inds){
      q = as.integer(infos$indicators$q[i])
      p = as.integer(infos$indicators$p[i])
      ind.lab = paste0("Y",q,".",p)
      # within
      loadW = mm_pars$sample[mm_pars$Level == "Within" & mm_pars$Type == "Loading"][i]
      loadW = ifelse(length(loadW) == 0, 1, loadW)
      sigmaW = mm_pars$sample[mm_pars$Level == "Within" & mm_pars$Type == "Measurement Error SD"][i]
      # between
      alpha = mm_pars$sample[mm_pars$Param == paste0("alpha_",q,".",p)]
      alpha = ifelse(length(alpha) == 0, 0, alpha)
      loadB = mm_pars$sample[mm_pars$Param == paste0("lambdaB_",q,".",p)]
      loadB = ifelse(length(loadB) == 0, 1, loadB)
      sigmaB = mm_pars$sample[mm_pars$Param == paste0("sigmaB_",q,".",p)]
      sigmaB = ifelse(length(sigmaB) == 0, 0, sigmaB)

      for(j in 1:N){
        # create WITHIN-PART:
        etaW = within[within$ID==j, paste0("Y",q)]
        if(sigmaW == 0){
          YW = loadW * etaW
        } else {
          YW = loadW * etaW + stats::rnorm(n = TP[j], mean = 0, sd = sigmaW)
        }

        # create BETWEEN-PART:
        if ( infos$p_is_wcen[i] == 1 ){
            YB = alpha + loadB * btw[j,infos$indicators$etaB_pos[i]] + ifelse(sigmaB == 0, 0, stats::rnorm(n = 1, mean = 0, sd = sigmaB))
          } else {
            YB = 0
          }
        # indicator
        within[within$ID==j,ind.lab] = YB + YW
      }
    }
    # remove latent process variables
    within = within[,!(colnames(within) %in% paste0("Y",1:infos$q))]
  }

  return(within)

}

# helper function to filter model object and add starting values for data generation
add_trues <- function(model, type = NULL, label = NULL, group = NULL, which = NULL, values, adjust_size = FALSE){

  select <- rep(TRUE, nrow(model))

  if(!is.null(type)){ select <- model$Type == type & select}
  if(!is.null(label)){ select <- model$Param_Label == label & select}
  if(!is.null(which)){ select <- which & select}
  if(!is.null(group)){ select <- model$group == group & select}


  if( adjust_size == TRUE ) { # fix length of values
    values = values[1:sum(select == TRUE)]
  }

  model$true.val[select] <- values

  return(model)
}


# reconstruct matrix of person-parameters of a fitted mlts.fit-object
get_person_par_mat <- function(
    fit,
    infos,
    iter){

  # construct matrix
  N = fit$standata$N
  g_id = fit$standata$g_id
  btw <- matrix(data = NA, nrow = N, ncol = infos$n_pars)

  # fill with samples
  ## random parameters
  n_rand = infos$n_random
  samples = rstan::extract(fit$stanfit, pars = paste0("b_free"))
  for(i in 1:n_rand){
    for(gg in 1:infos$G){
      btw[g_id==gg,infos$is_random[i]] <- samples$b_free[iter,g_id==gg,i]
    }
  }

  ## constants
  ### dynamic parameters
  n_fix = infos$n_fixed
  if(n_fix > 0){
    samples = rstan::extract(fit$stanfit, pars = paste0("b_fix"))
    for(i in 1:n_fix){
      for(gg in 1:infos$G){
      btw[g_id==gg,infos$is_fixed[i]] <- samples$b_fix[iter,gg,i]
      }
    }
  }
  ### innovation SDs
  n_innos_fix = infos$n_innos_fix
  if(n_innos_fix > 0){
    samples = rstan::extract(fit$stanfit, pars = paste0("sigma"))
    for(i in 1:n_innos_fix){
      for(gg in 1:infos$G){
        btw[g_id==gg,infos$innos_pos[i]] <- samples$sigma[iter,gg,i]
      }
    }
  }

  return(btw)
}


get_new_person_par_mat <- function(
    fit,
    infos,
    iter,
    gamma_pars,
    gamma_samples,
    sd_R_samples,
    bcorr_samples,
    re_pred_pars,
    re_pred_samples,
    W){


  # initial objects
  N = fit$standata$N
  g_id = fit$standata$g_id
  b_re_pred_mat <- list()
  bmu <- matrix(NA, nrow = N, ncol = infos$n_random)
  btw <- matrix(data = NA, nrow = N, ncol = infos$n_pars)
  n_rand <- infos$n_random
  b_free <- matrix(data = NA, nrow = N, ncol = n_rand)


  # loop over (potential groups)
  for(gg in 1:infos$G){
    b_re_pred_mat[[gg]] = matrix(NA, nrow = infos$n_cov, ncol = infos$n_random)
    for(j in 1:nrow(gamma_pars[gamma_pars$group==gg])){
      b_re_pred_mat[[gg]][1,] = gamma_samples$gammas[iter,gg,]
    }
    if(infos$n_cov>1){
      for(k in 1:infos$n_cov_bs){
        xx = infos$n_cov_mat[k,1]
        yy = infos$n_cov_mat[k,2]
        b_re_pred_mat[[gg]][xx,yy] = re_pred_samples[[re_pred_pars$Param_stan[re_pred_pars$group==gg][k]]][iter]
        }
      }

  # calculate population means
  bmu[g_id == gg,1:infos$n_random] = W[g_id == gg,] %*% b_re_pred_mat[[gg]]
  }

  # reconstruct var-cov-matrix
  for(p in 1:N){
    gg <- g_id[p]
    re_SDs = sd_R_samples$sd_R[iter,gg,]
    if(n_rand > 1){
      bcorr = as.matrix(bcorr_samples$bcorr[iter,gg,,])
      SIGMA = diag(re_SDs) %*% bcorr %*% diag(re_SDs)
      gammas = gamma_samples$gammas[iter,gg,]
      b_free[p,] = mvtnorm::rmvnorm(n = 1, mean = bmu[p,], sigma = SIGMA)

      } else {
      b_free[p,1] = bmu[p,1] + stats::rnorm(n = 1, mean = 0, sd = sd_R_samples$sd_R[iter,gg,1])
      }
    }
    for(i in 1:n_rand){
      btw[,infos$is_random[i]] <- b_free[,i]
    }

  ## constants
  ### dynamic parameters
  n_fix = infos$n_fixed
  if(n_fix > 0){
    samples = rstan::extract(fit$stanfit, pars = paste0("b_fix"))
    for(i in 1:n_fix){
      for(gg in 1:infos$G){
      btw[g_id == gg,infos$is_fixed[i]] <- samples$b_fix[iter,gg,i]
      }
    }
  }
  ### innovation SDs
  n_innos_fix = infos$n_innos_fix
  if(n_innos_fix > 0){
    samples = rstan::extract(fit$stanfit, pars = paste0("sigma"))
    for(i in 1:n_innos_fix){
      for(gg in 1:infos$G){
        btw[g_id == gg,infos$innos_pos[i]] <- samples$sigma[iter,gg,i]
      }
    }
  }

  return(btw)
}


