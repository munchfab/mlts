#' Title
#'
#' @param VARmodel data.frame. Output of VARmodel-Functions.
#' @param data data.frame. Data input.
#'
#' @return An object of class `data.frame`.
#' @export
#'
VARprepare <- function(VARmodel, data, ts.ind, covariates = NULL, outcomes = NULL,
                      center.covs = T, std.outcome = T
){


  # ignore most of data preprocessing steps already implemented in create_data ...


  # extract model inputs
  n_mus = sum(VARmodel$Param_Label == "Trait" & VARmodel$Type=="Fix effect")
  ## add check that n_mus == length(ts_ind)

  # subset of dynamic parameters (fixed effects)
  fix_pars = VARmodel[VARmodel$Type == "Fix effect",]
  fix_pars$no = 1:nrow(fix_pars)
  fix_pars_dyn = fix_pars[fix_pars$Param_Label == "Dynamic",]
  n_dyn.pars = nrow(fix_pars_dyn)
  if(n_mus < 10){  ### ACHTUNG AKTUELLE LÖSUNG FUNKTIONIERT NUR MIT D < 10
  fix_pars_dyn$Dout = substr(fix_pars_dyn$Param, start = 5, stop = 5)
  fix_pars_dyn$Dpred = substr(fix_pars_dyn$Param, start = 6, stop = 6)
  }
  # lagged relations between constructs
  N_pred = table(fix_pars_dyn$Dout) # number of lagged preds in each dimension
  D_pred = matrix(0, nrow = n_mus, ncol = n_mus, byrow = T)
  Dpos1 = c()
  Dpos2 = c()
  for(i in 1:n_mus){
    D_pred[i,1:N_pred[i]] = as.integer(fix_pars_dyn$Dpred[fix_pars_dyn$Dout == i])
    if(i == 1){
      Dpos1[i] = n_mus+1
      Dpos2[i] = n_mus+N_pred[i]
    } else {
      Dpos1[i] = Dpos2[i-1] +1
      Dpos2[i] = Dpos2[i-1] +N_pred[i]
    }
  }

  # which innovation variances as random effects?
  innos_rand = fix_pars[grepl(fix_pars$Param_Label, pattern="Variance"), "isRandom"]
  innos_pos = fix_pars[grepl(fix_pars$Param_Label, pattern="Variance"), "no"]
  n_innos_fix = n_mus - sum(innos_rand)
  innos_fix_pos = cumsum(1 - innos_rand)

  # number of innovation covariances to include
  n_inno_covs = nrow(fix_pars[grepl(fix_pars$Param_Label, pattern="Covariance"),])
  n_inno_cov_fix = n_inno_covs - sum(fix_pars[grepl(fix_pars$Param_Label, pattern="Covariance"),"isRandom"])
  inno_cov_pos = matrix(fix_pars[grepl(fix_pars$Param_Label, pattern="Covariance"), "no"],
                        nrow = 1, ncol = n_inno_covs)
  ##### ACHTUNG HIER AKTUELL NICHT MÖGLICH, DASS EINZELNE KOVARIANZEN ALS FIXED GESCHÄTZT WERDEN


  # parameters to pass to stan
  N = length(unique(data$ID))
  D = length(ts.ind)
  N_obs = nrow(data)
  N_obs_id = data.frame(table(data$ID))$Freq
  n_obs_cov = ifelse(n_inno_covs == 0, 0, N_obs - N) #### ACHTUNG!!!!!!

  n_pars = sum(VARmodel$Type == "Fix effect")
  n_random = sum(VARmodel$isRandom, na.rm = T)
  n_fixed = n_pars - n_random - n_innos_fix
  is_random = fix_pars$no[fix_pars$isRandom==1]
  is_fixed = matrix(fix_pars_dyn$no[fix_pars_dyn$isRandom==0], nrow = 1, ncol = n_fixed)



  ## handling of missing values
  n_miss = sum(is.na(data[,ts.ind]))
  n_miss_D = sapply(ts.ind, FUN = function(x){
    sum(is.na(data[,x]))
  })
  pos_miss_D = matrix(0, nrow = D, ncol = max(n_miss_D), byrow = T)
  for(i in 1:D){
    if(n_miss_D[i] > 0){
    pos_miss_D[i,1:n_miss_D[i]] = which(is.na(data[,ts.ind[i]]))
    }
  }

  ## replace missing values with -Inf
  data[,ts.ind][is.na(data[,ts.ind])] = -Inf
  y = t(data[,ts.ind])


  #### ADD BETWEEEN PART ....

  # Between-Level Covariates
  if(!is.null(covariates)){

    # check if number of covariates matches covariates in VARmodel
    covs = unique(model[model$Type == "Between-Level Regression","PRED"])
    n_cov = 1 + length(covs)  # add 1 for intercepts
    preds = model[model$Type == "Between-Level Regression",]
    random_pars = model[model$Type == "Fix effect" & model$isRandom==1,]
    n_cov_bs = nrow(preds)
    n_cov_mat = matrix(0, ncol = 2, nrow = n_cov_bs)
    for(i in 1:n_cov_bs){
      n_cov_mat[i,1] = which(covs == preds$PRED[i])+1
      n_cov_mat[i,2] = random_pars[random_pars$Param %in% preds$OUT[i], "par_number"]
    }
    W = matrix(NA, nrow = N, ncol = n_cov)
    W[,1] = 1
    for(i in 2:n_cov){
      for(p in 1:N){
        W[p,i] = unique(data[data$ID == p, covs[i-1]])
      }
      if(center.covs==T){
        W[,i] = W[,i] - mean(W[,i])
      }
    }

  } else {
    n_cov = 1
    n_cov_bs = n_random
    #### ACHTUNG NEEDS UPDATING
    n_cov_mat = matrix(1, ncol = n_cov+1, nrow = n_cov_bs)
    W = matrix(1, nrow = N, ncol = 1)
  }

  # Between-Level Outcomes
  if(!is.null(outcomes)){

    #### NEEDS UPDATING
    out_vars = unique(model[model$Type=="Outcome Prediction","OUT"])
    out_preds = unique(model[model$Type=="Outcome Prediction","PRED"])
    OUT = matrix(NA, nrow = length(out_vars), ncol = N)

    # out_pred_mat = matrix(NA, nrow = legnth(out_vars), ncol = )
    # if(length(out_vars)>0){
    #   for(p in 1:N){
    #     for(i in 1:length(out_vars)){
    #       W[i,p] = unique(data[data$ID == p, out_vars[i-1]])
    #       }
    #     }
    #   }
    # }

  }


  # get prior information from VARmodel
  if(n_random > 1){
    prior_LKJ = unique(VARmodel$prior_location[VARmodel$Type=="RE correlation"])
    prior_sigma_loc = VARmodel$prior_location[VARmodel$Type=="Random effect SD"]
    prior_sigma_scale = VARmodel$prior_scale[VARmodel$Type=="Random effect SD"]

    # add fix effect prior of constant innovation variance (as SD)
    if("Innovation Variance" %in% VARmodel$Param_Label){
      prior_sigma_loc = c(prior_sigma_loc,
                          VARmodel$prior_location[VARmodel$Param=="sigma_1"])
      prior_sigma_scale = c(prior_sigma_scale,
                            VARmodel$prior_scale[VARmodel$Param=="sigma_1"])
    }

    prior_gamma_loc = VARmodel$prior_location[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]
    prior_gamma_scale = VARmodel$prior_scale[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]

  } else {

    # NOT YET IMPLEMENTED
    prior_gamma_loc = VARmodel$prior_location[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]
    prior_gamma_scale = VARmodel$prior_scale[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]
    prior_sigma_loc = VARmodel$prior_location[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]
    prior_sigma_scale = VARmodel$prior_scale[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]
  }


  # combine all information
  standata = rstan::nlist(
    N, D, N_obs, N_obs_id, n_pars, n_random, n_fixed, is_random,
    is_fixed, y, n_miss, n_miss_D, pos_miss_D, innos_rand, innos_pos,
    n_innos_fix, innos_fix_pos, n_obs_cov, n_obs_cov, inno_cov_pos,
    n_cov, W, N_pred, D_pred, Dpos1, Dpos2, n_inno_covs, n_inno_cov_fix

    # priors
#    prior_LKJ, prior_sigma_loc, prior_sigma_scale, prior_gamma_loc, prior_gamma_scale
  )


  return(standata)

}
