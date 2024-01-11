#' Title
#'
#' @param VARmodel data.frame. Output of VARmodel-Functions.
#' @param data data.frame. Data input.
#'
#' @return An object of class `data.frame`.
#' @export
ARprepare <- function(VARmodel, data, ts.ind, covariates = NULL, outcomes = NULL,
                      center.covs = T, std.outcome = T
                      ){


  # ignore most of data preprocessing steps already implemented in create_data ...


  # evaluate model input to create data for stan
  N = length(unique(data$ID))
  N_obs = nrow(data)
  N_obs_id = data.frame(table(data$ID))$Freq
  n_pars = sum(VARmodel$Type == "Fix effect")
  n_random = sum(VARmodel$isRandom, na.rm = T)
  n_fixed = n_pars - n_random
  is_random = which(VARmodel$Type == "Fix effect" & VARmodel$isRandom==1)
  is_fixed = matrix(which(VARmodel$Type == "Fix effect" & VARmodel$isRandom==0),
                    ncol = n_fixed, nrow = 1)
  if(sum(VARmodel$Param_Label == "Log Innovation Variance")>0){
    innos_rand =  1
  } else {
    innos_rand = 0
  }


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


  y = unlist(data[,ts.ind])

  # handling of missing values
  n_miss = sum(data[,ts.ind] == -Inf)
  pos_miss = which(data[,ts.ind] == -Inf)


  # get prior information from VARmodel
  if(n_random > 1){
  prior_LKJ = unique(VARmodel$prior_location[VARmodel$Type=="Random effect correlation"])
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
    prior_LKJ = unique(VARmodel$prior_location[VARmodel$Type=="Random effect correlation"])
    prior_gamma_loc = VARmodel$prior_location[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]
    prior_gamma_scale = VARmodel$prior_scale[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]
    prior_sigma_loc = VARmodel$prior_location[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]
    prior_sigma_scale = VARmodel$prior_scale[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1]
  }


  # combine all information
  standata = rstan::nlist(
    N, N_obs, N_obs_id, n_pars, n_random, is_random,
    n_fixed, is_fixed, innos_rand, n_cov, n_cov_bs, n_cov_mat,
    W, y, n_miss, pos_miss,
    # priors
    prior_LKJ, prior_sigma_loc, prior_sigma_scale, prior_gamma_loc, prior_gamma_scale
  )


  return(standata)

}
