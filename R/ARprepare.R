#' Prepare list object of data to pass to `VARfit`-function for two-level AR(1) models
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


  # data entered here should be preprocessed prior (e.g., handling of missing values...)


  # data specific information:
  N = length(unique(data$num_id))                  # number of subjects
  N_obs = nrow(data)                           # total number of observations
  N_obs_id = data.frame(table(data$num_id))$Freq   # number of obs per subject

  ## replace missing values with -Inf
  data[,ts.ind][is.na(data[,ts.ind])] = -Inf
  y = unlist(data[,ts.ind])                    # store observations

  # handling of missing values
  n_miss = sum(data[,ts.ind] == -Inf)
  pos_miss = which(data[,ts.ind] == -Inf)


  # make use of VARmodelEval-function to extract information on VARmodel
  infos = VARmodelEval(VARmodel)


  n_pars = infos$n_pars       # no of dynamic parameters (including means, CRs, innovation variance)
  n_random = infos$n_random   # no of individual (random) effects
  n_fixed = infos$n_fixed     # no of fixed dynamic parameters (AR and CR effects)
  is_random = infos$is_random # which parameters (in order of Fix effect parameters in VARmodel) to model as random effect
  is_fixed = infos$is_fixed   # a matrix of n_fixed x 1, indicating which parameter to model as constant
  re_pars = infos$re_pars     # subset of VARmodel of random effect parameters
  innos_rand = infos$innos_rand # Are innovation variance(s) random?

  # between-level regression - REs as outcome
  RE.PREDS = infos$RE.PREDS
  n_cov = infos$n_cov  # 1 added for intercepts
  n_cov_bs = infos$n_cov_bs
  n_cov_mat = infos$n_cov_mat
  # build matrix of predictors
  W = matrix(NA, nrow = N, ncol = n_cov)
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
