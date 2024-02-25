#' Title
#'
#' @param VARmodel data.frame. Output of VARmodel-Functions.
#' @param data data.frame. Data input.
#'
#' @return An object of class `data.frame`.
#' @export
#'
VARprepare <- function(VARmodel, data, ts.ind, covariates = NULL, outcomes = NULL,
                       outcome.pred.btw = NULL, center.covs = T, std.outcome = T
){


  # data entered here should be preprocessed prior (e.g., handling of missing values...)

  # data specific information: -------------------------------------------------
  N = length(unique(data$num_id))         # number of subjects
  D = length(ts.ind)                      # number of time-varying constructs
  N_obs = nrow(data)                      # total number of observations (obs)
  N_obs_id = data.frame(table(data$num_id))$Freq # number of obs per subject

  ## handling of missing values
  n_miss = sum(is.na(data[,ts.ind]))           # overall number of NAs
  n_miss_D = sapply(ts.ind, FUN = function(x){ # store number of NAs per D
    sum(is.na(data[,x]))
  })
  n_miss_D = as.array(n_miss_D)
  # NA positions as matrix
  pos_miss_D = matrix(0, nrow = D, ncol = max(n_miss_D), byrow = T)
  for(i in 1:D){
    if(n_miss_D[i] > 0){
      pos_miss_D[i,1:n_miss_D[i]] = which(is.na(data[,ts.ind[i]]))
    }
  }

  ## replace missing values with -Inf
  data[,ts.ind][is.na(data[,ts.ind])] = -Inf
  y = t(data[,ts.ind])   # store observations

  # ----

  # model specific information: ------------------------------------------------
  infos = VARmodelEval(VARmodel)

  n_pars = infos$n_pars       # no of dynamic parameters (including means, CRs, innovation variance)
  n_random = infos$n_random   # no of individual (random) effects
  n_fixed = infos$n_fixed     # no of fixed dynamic parameters (AR and CR effects)
  is_random = as.array(infos$is_random) # which parameters (in order of Fix effect parameters in VARmodel) to model as random effect
  is_fixed = infos$is_fixed   # a matrix of n_fixed x 1, indicating which parameter to model as constant
  re_pars = infos$re_pars     # subset of VARmodel of random effect parameters

  # innovations
  innos_rand = as.array(infos$innos_rand) # Are innovation variance(s) random?
  n_innos_fix = infos$n_innos_fix
  innos_fix_pos = as.array(infos$innos_fix_pos)
  innos_pos = as.array(infos$innos_pos)

  # innovation covariances
  n_inno_covs = infos$n_inno_covs
  n_inno_cors = infos$n_inno_cors
  n_inno_cov_fix = infos$n_inno_cov_fix
  inno_cov_pos = infos$inno_cov_pos

  # NEEDS TO BE UPDATED ----------------------------
  n_obs_cov = ifelse(n_inno_covs == 0, 0, N_obs - N)
  # -----------------------------------------


  # dynamic model specification by dimension (D)
  maxLag = infos$maxLag
  N_pred = infos$N_pred
  D_pred = infos$D_pred
  Lag_pred = infos$Lag_pred
  Dpos1 = as.array(infos$Dpos1)
  Dpos2 = as.array(infos$Dpos2)
  # covariate(s) as predictor(s) of random effects
  RE.PREDS = infos$RE.PREDS
  n_cov = infos$n_cov
  n_cov_bs = infos$n_cov_bs
  n_cov_mat = infos$n_cov_mat
  W = matrix(NA, nrow = N, ncol = n_cov) # build matrix of predictors
  W[,1] = 1
  if(n_cov > 1){
    for(i in 2:n_cov){
      cov.use = which(covariates == unique(RE.PREDS$re_preds[RE.PREDS$pred_no == i-1]))
      for(p in 1:N){
        W[p,i] = unique(data[data$num_id == p, names(covariates)[cov.use]])
      }
      if(center.covs==T){
        W[,i] = W[,i] - mean(W[,i])
      }
    }
  }
  # outcome prediction - REs as predictors
  n_out = infos$n_out
  out_vars = infos$out_var
  out_vars_data = unlist(lapply(out_vars, function(x){
    names(outcomes)[which(outcomes == x)]
  }))
  n_out_bs = infos$n_out_bs
  n_out_b_pos = infos$n_out_b_pos
  n_out_bs_max = infos$n_out_bs_max
  n_out_bs_sum = infos$n_out_bs_sum
  # outcome data
  out = matrix(nrow = n_out, ncol = N, data = NA)
  if(n_out>0){
    for(i in 1:n_out){
      for(p in 1:N){
        out[i,p] = unique(data[data$num_id==p,out_vars_data[i]])
      }
    }
  }

  n_z_vars = infos$n_z_vars
  n_z_vars_data = unlist(lapply(n_z_vars, function(x){
    names(outcome.pred.btw)[which(outcome.pred.btw == x)]
  }))
  n_z = infos$n_z
  Z = matrix(NA, nrow = N, ncol = n_z)
  if(n_z > 0){
    for(i in 1:n_z){
      for(p in 1:N){
        Z[p,i] = unique(data[data$num_id==p, n_z_vars_data[i]])
      }
    }
  }

  # priors
  prior_gamma = infos$prior_gamma
  prior_sd_R = infos$prior_sd_R
  prior_LKJ = infos$prior_LKJ
  prior_sigma = infos$prior_sigma
  prior_b_re_pred = infos$prior_b_re_pred
  prior_b_out = infos$prior_b_out
  prior_alpha_out = infos$prior_alpha_out
  prior_sigma_out = infos$prior_sigma_out



  # combine all information
  standata = rstan::nlist(
    # data specifications
    N, D,maxLag, N_obs, N_obs_id, y, n_miss, n_miss_D, pos_miss_D,
    # model specifications
    n_pars, n_random, n_fixed, is_random, is_fixed,
    N_pred, D_pred, Lag_pred, Dpos1, Dpos2,
    ## innovations
    innos_rand, innos_pos, n_innos_fix,
    ## innovation covariances
    innos_fix_pos, n_obs_cov, n_obs_cov, inno_cov_pos,
    n_inno_covs,n_inno_cors, n_inno_cov_fix,

    n_cov, n_cov_bs, n_cov_mat, W,
    n_out, n_out_bs, n_out_bs_max, n_out_bs_sum, n_out_b_pos, out, n_z, Z,
    # priors
    prior_LKJ, prior_gamma, prior_sd_R, prior_sigma, prior_b_re_pred,
    prior_b_out, prior_alpha_out, prior_sigma_out
  )

  # for latent models:
  if(infos$isLatent == TRUE){
    standata$D = infos$q
    standata$n_p = nrow(infos$indicators)
    standata$D_perP = as.array(as.integer(infos$indicators$q))

    # alternative missing data indexing
    n_miss = sum(data[,ts.ind] == -Inf)           # overall number of NAs
    n_miss_p = as.array(sapply(ts.ind, FUN = function(x){sum(data[,x] == -Inf)}))
    pos_miss_p = matrix(0, nrow = standata$n_p, ncol = max(n_miss_p), byrow = T)
    for(i in 1:standata$n_p){
      if(n_miss_p[i] > 0){
        pos_miss_p[i,1:n_miss_p[i]] = which(data[,ts.ind[i]] == -Inf)
      }
    }
    standata$n_miss <- n_miss
    standata$n_miss_p <- n_miss_p
    standata$pos_miss_p <- pos_miss_p

    standata$n_loadBfree <- infos$n_loadBfree
    standata$n_loadWfree <- infos$n_loadWfree
    standata$n_alphafree <- infos$n_alphafree
    standata$n_sigmaBfree <- infos$n_sigmaBfree
    standata$n_sigmaWfree <- infos$n_sigmaWfree
    standata$pos_loadBfree <- as.array(infos$pos_loadBfree)
    standata$pos_loadWfree <- as.array(infos$pos_loadWfree)
    standata$pos_alphafree <- as.array(infos$pos_alphafree)
    standata$pos_sigmaBfree <- as.array(infos$pos_sigmaBfree)
    standata$pos_sigmaWfree <- as.array(infos$pos_sigmaWfree)

    standata$n_YB_free <- infos$n_YB_free
    standata$YB_free_pos <- as.array(infos$YB_free_pos)
    standata$mu_etaB_pos <- as.array(infos$mu_etaB_pos)
    standata$mu_is_etaB <- as.array(infos$mu_is_etaB)
  }


  return(standata)

}
