#' Title
#'
#' @param VARmodel data.frame. Output of VARmodel-Functions.
#' @param data data.frame. Data input.
#'
#' @return An object of class `data.frame`.
#' @export
#'
ARprepare <- function(VARmodel, data, ts.ind, covariates = NULL, outcomes = NULL,
                      center.covs = T, std.outcome = T, outcome.pred.btw = NULL
                      ){

  # ignore most of data preprocessing steps already implemented in create_data ...

  # evaluate model input to create data for stan
  infos = VARmodelEval(VARmodel)

  N = length(unique(data$ID))
  N_obs = nrow(data)
  N_obs_id = data.frame(table(data$ID))$Freq
  n_pars = sum(VARmodel$Type == "Fix effect")
  n_random = sum(VARmodel$isRandom, na.rm = T)
  n_fixed = n_pars - n_random
  is_random = which(VARmodel$Type == "Fix effect" & VARmodel$isRandom==1)
  is_fixed = matrix(which(VARmodel$Type == "Fix effect" & VARmodel$isRandom==0),
                    ncol = n_fixed, nrow = 1)
  re_pars = VARmodel[VARmodel$Type=="Fix effect" & VARmodel$isRandom==1,]
  re_pars$par_no = 1:nrow(re_pars)
  if(sum(VARmodel$Param_Label == "Log Innovation Variance")>0){
    innos_rand =  1
  } else {
    innos_rand = 0
  }

  # evaluate model inputs for between-level regressions
  ## REs as outcomes
  RE.PREDS = VARmodel[VARmodel$Type == "RE prediction",]

  if(nrow(RE.PREDS)>0  | !is.null(covariates)){
    # which REs to regress on
    RE.PREDS$re_as_dv = substring(
      RE.PREDS$Param, 3,regexpr(RE.PREDS$Param,pattern = ".ON.", fixed = T)-1)
    RE.PREDS$re_preds = substring(
      RE.PREDS$Param, regexpr(RE.PREDS$Param,pattern = ".ON.", fixed = T)+4)
    RE.PREDS$re_preds_model = sapply(RE.PREDS$re_preds, function(x){
      names(covariates)[which(covariates == x)]})
    RE.PREDS$re_pred_no = sapply(RE.PREDS$re_preds, function(x){
      which(covariates == x)})
    RE.PREDS$re_no = sapply(RE.PREDS$re_as_dv, function(x){
      re_pars$par_no[which(re_pars$Param == x)]})
    RE.PREDS$re_pred_b_no = 1:nrow(RE.PREDS)

    re_preds_unique = unique(RE.PREDS$re_preds)
    # check if covariates in VARmodel match covariates in the data
    if(sum(re_preds_unique %in% covariates) != length(re_preds_unique)){
      warning("Between-level covariate names in VARmodel and
               covariate argument do not match!")
    }

    n_cov = 1 + length(covariates)                        # add 1 for intercepts
    n_cov_bs = nrow(RE.PREDS)
    n_cov_mat = matrix(unlist(RE.PREDS[,c("re_pred_no", "re_no")]),
                       ncol = 2, nrow = n_cov_bs)
    n_cov_mat[,1] = n_cov_mat[,1] + 1    # shift by 1 for intercepts
    W = matrix(NA, nrow = N, ncol = n_cov)
    W[,1] = 1
    for(i in 2:n_cov){
      for(p in 1:N){
        W[p,i] = unique(data[data$ID == p, names(covariates)[i-1]])
      }
      if(center.covs==T){
        W[,i] = W[,i] - mean(W[,i])
      }
    }
  } else { # covariates
    n_cov = 1
    n_cov_bs = 0    # use a random placeholder as it will be overwritten when n_cov = 1
    n_cov_mat = matrix(1, ncol = n_cov+1, nrow = n_cov_bs)
    W = matrix(1, nrow = N, ncol = 1)
  }


  # evaluate model inputs for outcome prediction
  ## REs as outcomes
  OUT = VARmodel[VARmodel$Type == "Outcome prediction" & VARmodel$Param_Label == "regression weight",]
  if(nrow(OUT)>0 | !is.null(outcomes)){
    OUT$out = substring(OUT$Param, 3,regexpr(OUT$Param,pattern = ".ON.", fixed = T)-1)
    OUT$out_pred = substring(OUT$Param, regexpr(OUT$Param,pattern = ".ON.", fixed = T)+4)
    OUT$out_model = sapply(OUT$out, function(x){
      names(outcomes)[which(outcomes == x)]})

    # extract info whether additional btw pars are included
    n_z_vars = unique(OUT$out_pred[!(OUT$out_pred%in% re_pars$Param)])
    n_z = length(n_z_vars)
    Z = matrix(NA, nrow = N, ncol = n_z)
    for(i in 1:n_z){
      for(p in 1:N){
        Z[p,i] = unique(data[data$ID==p, names(outcome.pred.btw[i])])
      }
    }
    OUT$out_pred_no = sapply(OUT$out_pred, function(x){
      ifelse(x %in% re_pars$Param, re_pars$par_no[which(re_pars$Param==x)],
             which(outcome.pred.btw == x) + n_random)})

    out_vars = unique(OUT$out)
    n_out = length(out_vars)
    n_out_bs = matrix(ncol = 1,nrow = n_out)
    for(i in 1:n_out){
      n_out_bs[i,1] = sum(OUT$out == out_vars[i])
    }
    n_out_b_pos = matrix(0, ncol = max(n_out_bs[,1]), nrow = n_out)
    for(i in 1:n_out){
      n_out_b_pos[i,1:n_out_bs[i,1]] = OUT$out_pred_no[OUT$out==out_vars[i]]
    }
    n_out_bs_sum = sum(n_out_bs[,1])
    n_out_bs_max = max(n_out_bs[,1])

    # outcome data
    out = matrix(nrow = n_out, ncol = N, data = NA)
    for(i in 1:n_out){
      for(pp in 1:N){
        out[i,pp] = unique(data[data$ID==pp,names(outcomes[which(outcomes == out_vars[i])])])
      }
    }

  } else{
    n_out = 0
    n_out_bs = matrix(ncol = 1, nrow = n_out)
    n_out_bs_sum = 0
    n_out_bs_max = 0
    n_out_b_pos = matrix(nrow = n_out, ncol = 0)
    out = matrix(nrow = 0, ncol = N)
    n_z = 0
    Z = matrix(nrow = N, ncol = 0)
  }

  y = unlist(data[,ts.ind])

  # handling of missing values
  n_miss = sum(data[,ts.ind] == -Inf)
  pos_miss = which(data[,ts.ind] == -Inf)


  # get prior information from VARmodel
  if(n_random > 1){

    # fixed effects (intercepts)
    prior_gamma = matrix(nrow = n_random, ncol = 2)
    prior_gamma[,1] = VARmodel$prior_location[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]
    prior_gamma[,2] = VARmodel$prior_scale[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]

    # random effect SDs
    prior_sd_R = matrix(nrow = n_random, ncol = 2)
    prior_sd_R[,1] = VARmodel$prior_location[VARmodel$Type=="Random effect SD"]
    prior_sd_R[,2] = VARmodel$prior_scale[VARmodel$Type=="Random effect SD"]

    # random effect correlations
    prior_LKJ = unique(VARmodel$prior_location[VARmodel$Type=="RE correlation"])

    # add fix effect prior of constant innovation variance (as SD)
    prior_sigma = matrix(ncol = 2, nrow = sum(infos$n_innos_fix))
    if(innos_rand == 0){
      prior_sigma[1,1] = VARmodel$prior_location[VARmodel$Param=="sigma_1"]
      prior_sigma[1,2] = VARmodel$prior_scale[VARmodel$Param=="sigma_1"]
    }

    # RE prediction
    prior_b_re_pred = matrix(ncol = 2, nrow = n_cov_bs)
    if(!is.null(covariates)){
      prior_b_re_pred[,1] = RE.PREDS$prior_location
      prior_b_re_pred[,2] = RE.PREDS$prior_scale
    }

    # outcome prediction
    prior_b_out = matrix(ncol = 2, nrow = n_out_bs_sum)
    prior_alpha_out = matrix(ncol = 2, nrow = n_out)
    prior_sigma_out = matrix(ncol = 2, nrow = n_out)

    if(!is.null(outcomes)){
      infos$OUT = infos$OUT[order(infos$OUT$out_var_no, infos$OUT$Pred_no),]
      prior_b_out[,1] = infos$OUT$prior_location
      prior_b_out[,2] = infos$OUT$prior_scale
      for(i in 1:n_out){
        prior_alpha_out[i,1] = VARmodel$prior_location[VARmodel$Param %in% c(paste0("alpha_",infos$out_var[i]))]
        prior_alpha_out[i,2] = VARmodel$prior_scale[VARmodel$Param %in% c(paste0("alpha_",infos$out_var[i]))]
        prior_sigma_out[i,1] = VARmodel$prior_location[VARmodel$Param %in% c(paste0("sigma_",infos$out_var[i]))]
        prior_sigma_out[i,2] = VARmodel$prior_scale[VARmodel$Param %in% c(paste0("sigma_",infos$out_var[i]))]
      }
    }

  } else {

    # NOT YET IMPLEMENTED
    prior_LKJ = unique(VARmodel$prior_location[VARmodel$Type=="RE correlation"])
    prior_gamma = VARmodel$prior_location[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]
    prior_gamma = VARmodel$prior_scale[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]
    prior_sd_R = VARmodel$prior_location[VARmodel$Type=="Random effect SD"]
    prior_sd_R = VARmodel$prior_scale[VARmodel$Type=="Random effect SD"]
    prior_sigma = VARmodel$prior_location[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]
    prior_sigma = VARmodel$prior_scale[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]
    prior_b_re_pred = matrix(RE.PREDS$prior_location, nrow = 1, ncol = n_cov_bs)
    prior_b_re_pred = matrix(RE.PREDS$prior_scale, nrow = 1, ncol = n_cov_bs)
  }


  # combine all information
  standata = rstan::nlist(
    N, N_obs, N_obs_id, n_pars, n_random, is_random,
    n_fixed, is_fixed, innos_rand, n_cov, n_cov_bs, n_cov_mat,
    W, y, n_miss, pos_miss,
    n_out, n_out_bs, n_out_bs_max, n_out_bs_sum, n_out_b_pos, out, n_z, Z,
    # priors
    prior_LKJ, prior_gamma, prior_sd_R, prior_sigma, prior_b_re_pred,
    prior_b_out, prior_alpha_out, prior_sigma_out
  )


  return(standata)

}
