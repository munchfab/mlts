#' Simulate data from mlts model
#'
#' @details
#' A function to generate data from an output of \code{\link[mlts]{mlts_model}}.
#'
#' @param model `data.frame`. Output of \code{\link[mlts]{mlts_model}}.
#' @param default logical. If set to `TRUE`, default prior specifications are
#' added.
#' @param N integer Number of observational units.
#' @param TP integer. Number of measurements per observational unit.
#' @param burn.in integer. Length of ‘burn-in’ period.
#' @param seed integer. Seed used for data generation.
#' @param seed.true integer. Separate seed used for sampling of true
#' population parameters values from plausible ranges for stationary time series.
#' @param btw.var.sds named numeric vector. Provide standard deviation(s) for all exogenous
#' between-level variable(s) specified in `model`, e.g. (`btw.var.sds = c("covariate1" = 1)`,
#' to set the SD of the variable "covariate1" to 1). Mean values of the respective
#' variable(s) will be set to 0 per default.
#' @return An object of class \code{"mlts_simdata"}.
#' The object is a list containing the following components:
#' \item{model}{the model object passed to `mlts_sim` with true parameter values used
#' in the data generation added in the column `true.val`}
#' \item{data}{a long format `data.frame` of the generated time series data}
#' \item{RE.pars}{a `matrix` of cluster-specific true values used in the data generation}
#' @export
#'
#' @examples
#' \donttest{
#' # build a simple vector-autoregressive mlts model with two time-series variables
#' var_model <- mlts_model(q = 2)
#'
#' # simulate data from this model with default true values
#' # (true values are randomly drawn from normal distribution)
#' var_data <- mlts_sim(
#'   model = var_model,
#'   N = 50, TP = 30, # number of units and number of measurements per unit
#'   default = TRUE # use default parameter values
#' )
#'
#' # the data set is stored in .$data
#' head(var_data$data)
#'
#' # individual parameter values are stored in .$RE.pars
#' head(var_data$RE.pars)
#'
#' # if the mltssim-object is used in mlts_fit(), true values
#' # are added to the fitted object
#' fit <- mlts_fit(
#'   model = var_model,
#'   data = var_data,
#'   id = "ID", ts = c("Y1", "Y2"), time = "time"
#' )
#'
#' # inspect model with true values
#' head(fit$pop.pars.summary)
#' }
#'

mlts_sim <- function(model, default = FALSE, N, TP, burn.in = 50, seed = NULL,
                     seed.true = 1, btw.var.sds = NULL){

  set.seed(seed.true)

  ### start with a check that sds are entered for all between variables


  # simulate data based on model

  ### write separate function for specification of true parameter values
  ### for now use fixed values


  # use helper function to read out information on model
  infos = mlts_model_eval(model)


  # use some default settings for parameter values
  if(default==TRUE){
    model$true.val = NA

    # MEASUREMENT MODEL PARAMETERS =============================================
    n_p = sum(infos$p)
    Model = "Measurement"
    model$true.val[model$Model == Model & model$Type == "Item intercepts" & model$Constraint == "= 0"] = 0
    model$true.val[model$Model == Model & model$Type == "Item intercepts" & model$Constraint == "free"] =
      sample(x = seq(0.5, 2, by = 0.5), size = infos$n_alphafree, replace = TRUE)
    model$true.val[model$Level == "Within" & model$Type == "Loading" & model$Constraint == "= 1"] = 1
    model$true.val[model$Level == "Within" & model$Type == "Loading" & model$Constraint == "free"] =
      sample(x = seq(0.7, 0.9, by = 0.05), size = infos$n_loadWfree, replace = TRUE)
    model$true.val[model$Level == "Between" & model$Type == "Loading" & model$Constraint == "= 1"] = 1
    model$true.val[model$Level == "Between" & model$Type == "Loading" & model$Constraint == "free"] =
      sample(x = seq(0.7, 0.9, by = 0.05), size = infos$n_loadBfree, replace = TRUE)
    model$true.val[model$Level == "Within" & model$Type == "Measurement Error SD" & model$Constraint == "= 0"] = 0
    model$true.val[model$Level == "Within" & model$Type == "Measurement Error SD" & model$Constraint == "free"] =
      sample(x = seq(0.15, 0.3, by = 0.05), size = infos$n_sigmaWfree, replace = TRUE)
    model$true.val[model$Level == "Between" & model$Type == "Measurement Error SD" & model$Constraint == "= 0"] = 0
    model$true.val[model$Level == "Between" & model$Type == "Measurement Error SD" & model$Constraint == "free"] =
      sample(x = seq(0.15, 0.2, by = 0.05), size = infos$n_sigmaBfree, replace = TRUE)

    # Fixed effects ========
    model.type = "Fixed effect"
    ## Mus
    n_traits = length(model$true.val[model$Type==model.type & startsWith(model$Param_Label, "Trait")])
    model$true.val[model$Type==model.type & startsWith(model$Param_Label, "Trait")] = sample(
      x = seq(from= 0, to = 1, by = 0.1), size = n_traits)
    ## Phis
    # dynamic parameters separately for each lag
    phis = infos$fix_pars_dyn
    phis$true.val = NA
    ## start with lag 1 ARs
    n_sample = nrow(phis[phis$Lag == 1 & phis$isAR == 1,])
    phis$true.val[phis$Lag == 1 & phis$isAR == 1] = sample(x = seq(from=.15,to=0.3,by=0.05), replace = TRUE, size = n_sample)
    n_sample = nrow(phis[phis$Lag == 1 & phis$isAR == 0,])
    phis$true.val[phis$Lag == 1 & phis$isAR == 0] = sample(x = seq(from=-0.2,to=0.1,by=0.05),replace = TRUE, size = n_sample)
    ## add higher order effects
    if(infos$maxLag>1){
      for(i in 1:nrow(phis)){
        if(phis$Lag[i] > 1){
          phi.lag1 = phis$true.val[phis$Param == paste0("phi(1)_",phis$Dout[i],phis$Dpred[i])]
          phis$true.val[i] = phi.lag1/phis$Lag[i]
        }
      }
    }
    model$true.val[model$Type==model.type & model$Param_Label=="Dynamic"] = round(phis$true.val,3)


    ## log innovation variances
    model$true.val[model$Type==model.type & model$Param_Label=="Log Innovation Variance"] = -0.3
    ## Fixed innovation variance
    model$true.val[model$Type==model.type & model$Param_Label=="Innovation Variance"] = 0.75

    ## Log innovation covariance
    model$true.val[model$Type==model.type & model$Param_Label=="Log Innovation Covariance"] = -0.3
    model$true.val[model$Type==model.type & model$Param_Label=="Innovation correlation"] = -0.15

    # RANDOM EFFECT SDs ==========
    model.type = "Random effect SD"
    ## Mus
    model$true.val[model$Type==model.type & startsWith(model$Param_Label, "Trait")] = sample(
      x = seq(from= 0.7, to = 1.2, by = 0.1), size = n_traits, replace = TRUE)
    ## Phis
    model$true.val[model$Type==model.type & model$Param_Label=="Dynamic"] = 0.15
    ## log innovation variances
    model$true.val[model$Type==model.type & model$Param_Label=="Log Innovation Variance"] = 0.25
    ## log innovation covaraince(s)
    model$true.val[model$Type==model.type & model$Param_Label=="Log Innovation Covariance"] = 0.25


    # RANDOM EFFECT CORRELATIONS ============
    ## set all to zero for now
    model.type = "RE correlation"
    model$true.val[model$Type == model.type] = sample(
      c(round(seq(from = -0.2, to = 0.2, by = 0.025),3)), replace = TRUE,
      size = sum(model$Type == model.type))

    # RE as OUTCOME =========================
    model.type = "RE prediction"
    model$true.val[model$Type == model.type] = sample(
      c(round(seq(from = -0.2, to = 0.2, by = 0.05),3)), replace = TRUE,
      size = sum(model$Type == model.type))

    # OUTCOME PREDICTION ===================
    model.type = "Outcome prediction"
    model$true.val[model$Type == model.type] = sample(
      c(round(seq(from = -0.3, to = 0.3, by = 0.1),3)), replace = TRUE,
      size = sum(model$Type == model.type))
    # scale true values for AR and CL as predictor
    model$true.val[model$Type == model.type & grepl(pattern = "phi",model$Param)] <-
      model$true.val[model$Type == model.type & grepl(pattern = "phi",model$Param)] * 5

    model$true.val[model$Type == model.type & model$Param_Label == "intercept"] = 0
    model$true.val[model$Type == model.type & model$Param_Label == "Residual SD"] = 0.5

  } else if ( is.null(model$true.val)) {
    stop("No true parameter values provided in model$true.val. Set default = TRUE to run data generation with random true parameter values.",
         "Alternatively, user-specified values for each parameter can be specified in an additional column `true.val` in model.")
  }


  # set seed for data generation
  if(!is.null(seed)){
    set.seed(seed)
  }



  # run again after adding true parameter values
  infos = mlts_model_eval(model)


  # start generating between-level model =======================================

  # FIXED EFFECTS
  gammas = model$true.val[model$Type=="Fixed effect" & model$isRandom==1]

  # BETWEEN-LEVEL
  # sample covariates and get expected values of individual parameters
  bmu = matrix(data = NA, nrow = N, ncol = infos$n_random)
  W = matrix(data = NA, nrow = N, ncol = infos$n_cov)
  cov_name = c()
  W[,1] = 1 # intercept
  if(infos$n_cov>1){
    for(i in 2:infos$n_cov){
      cov_name[i-1] = unique(infos$RE.PREDS$re_preds[infos$RE.PREDS$pred_no == i-1])
      # name btw.var.sds vector
      # names(btw.var.sds) <- unique(infos$RE.PREDS$re_preds[infos$RE.PREDS$pred_no == i-1])
      W[1:N,i] = stats::rnorm(n = N, mean = 0, btw.var.sds[names(btw.var.sds) == cov_name[i-1]])
    }
  }
  colnames(W) <- c("Intercept", cov_name)

  for(i in 1:infos$n_random){
    # get expected individual parameters
    pred_use = infos$RE.PREDS[infos$RE.PREDS$re_no ==i,]
    if(nrow(pred_use)>0){
      bmu[,i] = W[,c(1,pred_use$pred_no+1)] %*% c(gammas[i],pred_use$true.val)
    } else {
      bmu[,i] = gammas[i]
    }
  }


  # calculate covariances from correlations
  n_random = infos$n_random
  rand.pars = infos$re_pars$Param

  # variance covariance matrix of random effects
  if(n_random == 1){
    cov_mat = model$true.val[model$Type=="Random effect SD"]
  } else {
  cov_mat = diag(model$true.val[model$Type=="Random effect SD"]^2)
  for(i in 1:n_random){
    for(j in 1:n_random){
      if(i < j){
        r = model$true.val[model$Param == paste0("r_",rand.pars[i],".", rand.pars[j])]
        cov_mat[i,j] = cov_mat[j,i] <- r * sqrt(cov_mat[i,i]) * sqrt(cov_mat[j,j])
        }
      }
    }
  }


  #### sample random effects from multivariate normal distribution and add to bmus
  btw_random = matrix(NA, nrow = N, ncol = infos$n_random)
  if(n_random == 1){
    btw_random = bmu + stats::rnorm(n = N, mean = 0, sd = cov_mat)
  } else {
    btw_random = bmu + mvtnorm::rmvnorm(n = N, mean = rep(0, infos$n_random), sigma = cov_mat)
  }
  colnames(btw_random) = infos$fix_pars$Param[infos$is_random]

  # check for AR parameters with absolute values below "1"
  posAR = infos$fix_pars[infos$fix_pars$isRandom==1 & infos$fix_pars$isAR==1, "no"]

  if(sum(abs(btw_random[,posAR]) >= 1)){
    stop("Absolute individual AR effects greater than 1 were sampled. Consider
            setting the true values of the fixed effect or the random effect SD to a lower value.")
  }

  # now combine fixed effects and random effects
  btw = matrix(NA, nrow = N, infos$n_pars)
  btw[,infos$is_random] = btw_random               # first add random pars
  if(infos$n_fixed>0){
    btw[,infos$is_fixed[1,]] = rep(model$true.val[model$Type=="Fixed effect"][infos$is_fixed[1,]], each=N)
  }
  if(infos$n_innos_fix>0){
    for(i in infos$innos_fix_pos)
      btw[,infos$innos_pos[i]] = rep(model$true.val[model$Type=="Fixed effect" & model$Param_Label == "Innovation Variance"][i],times=N)
  }

  #### WITHIN-LEVEL PROCESS ====================================================
  # prepare data frame
  NT = TP + burn.in
  within = data.frame(
    "ID" = rep(1:N, each = TP),
    "time" = rep(1:TP, times = N)
  )
  q = infos$q   # number of constructs
  y_cols = paste0("Y",1:q)   # prepare columns
  within[, y_cols] = NA

  # get positions to fill transition matrix
  dyn = infos$fix_pars_dyn
  dyn$tran_pos = NA
  dyn$tran_pos = as.integer(dyn$Dpred) + q*(dyn$Lag-1)

  for(i in 1:N){
    # build person-specific transition matrix
    transition = matrix(nrow = q, ncol = q*infos$maxLag, data = 0)
    y = matrix(nrow = NT, ncol = q)

    # build person-specific prediction error matrix
    innoVars.i = btw[i,infos$innos_pos]
    for(d in 1:q){ # loop over number of constructs

      transition[d,as.integer(dyn$tran_pos[dyn$Dout==d])] = btw[i,infos$Dpos1[d]:infos$Dpos2[d]]
      if(infos$innos_rand[d] == 1){
        innoVars.i[d] = exp(innoVars.i[d]) # retransform log innovations
      } else{
        innoVars.i[d] = innoVars.i[d]^2
      }
    }

    if(q == 1){
      inno_var_mat = matrix(data = innoVars.i, nrow = 1, ncol = 1)
    # } else if(q == 2 & infos$n_inno_covs == 1){
    #   inno_var_mat = diag(innoVars.i)
    #   inno_var_mat[1,2] <- inno_var_mat[2,1] <- exp(btw[i, infos$inno_cov_pos])
    } else {
      inno_var_mat = diag(innoVars.i)

      if(infos$n_inno_cors > 0){
        for(xx in 1:q){
          for(yy in 1:q){
            if(xx < yy){
              cor = model$true.val[model$Param == paste0("r.zeta_",xx,yy)]
              cov = cor * sqrt(inno_var_mat[xx,xx]) * sqrt(inno_var_mat[yy,yy])
              inno_var_mat[xx,yy] <- inno_var_mat[yy,xx] <- cov
            }
          }
        }
      }
    }

    for(t in 1:NT){
      # use means as starting values
      if(t <= infos$maxLag){
        init = mvtnorm::rmvnorm(n = 1, mean = rep(0, q), sigma = inno_var_mat)
        y[t,] = init
      } else {
        y_lag = c()
        y_lag = as.vector(y[t-1,1:q])
        if(infos$maxLag>1){
          for(ll in 2:infos$maxLag){
            y_lag = c(y_lag, as.vector(y[t-ll,1:q]))
          }
        }

        if(q > 1 | infos$maxLag > 1){
          y[t,1:q] = y_lag %*% t(transition) # get expected values
        } else {
          y[t,] = y_lag * transition
        }
        y[t,] = y[t,] + mvtnorm::rmvnorm(n = 1, mean = rep(0, q), sigma = inno_var_mat) # add error


        if(infos$q == 2 & infos$n_inno_covs == 1){
          inno_t = stats::rnorm(n = 1, mean = 0, sd = sqrt(exp(btw[i,infos$inno_cov_pos])))
          y[t,] = y[t,] + infos$inno_cov_load * inno_t
        }
      }
    }
    within[within$ID==i, y_cols] = y[(burn.in+1) : (burn.in+TP),]

    # add trait scores (for manifest indicators)
    if(infos$isLatent == FALSE){
      for(j in 1:infos$q){
        within[within$ID==i,y_cols[j]] = btw[i,j] + within[within$ID==i,y_cols[j]]
      }
    }
  }

  # # remove burn.in
  # within$time = within$time - burn.in
  # within = within[within$time>0,]
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
      loadW = model$true.val[model$Level == "Within" & model$Type == "Loading"][i]
      loadW = ifelse(length(loadW) == 0, 1, loadW)
      sigmaW = model$true.val[model$Level == "Within" & model$Type == "Measurement Error SD"][i]
      # between
      alpha = model$true.val[model$Param == paste0("alpha_",q,".",p)]
      alpha = ifelse(length(alpha) == 0, 0, alpha)
      loadB = model$true.val[model$Param == paste0("lambdaB_",q,".",p)]
      loadB = ifelse(length(loadB) == 0, 1, loadB)
      sigmaB = model$true.val[model$Param == paste0("sigmaB_",q,".",p)]
      sigmaB = ifelse(length(sigmaB) == 0, 0, sigmaB)

      for(j in 1:N){
        # create WITHIN-PART:
        etaW = within[within$ID==j, paste0("Y",q)]
        if(sigmaW == 0){
          YW = loadW * etaW
        } else {
          YW = loadW * etaW + stats::rnorm(n = TP, mean = 0, sd = sigmaW)
        }

        # create BETWEEN-PART:
        YB = alpha + loadB * btw[j,infos$indicators$etaB_pos[i]] + ifelse(sigmaB == 0, 0, stats::rnorm(n = 1, mean = 0, sd = sigmaB))

        # indicator
        within[within$ID==j,ind.lab] = YB + YW
      }
    }
    # remove latent process variables
    within = within[,!(colnames(within) %in% paste0("Y",1:infos$q))]
  }

  # CREATE OUTCOMES ==========================================================
  outs = matrix(NA, ncol = infos$n_out, nrow = N)
  if(infos$n_out > 0){
    if(infos$n_z > 0){
      btw.Z = matrix(nrow = N, ncol = (infos$n_random + infos$n_z))
      btw.Z[,1:infos$n_random] = btw_random
      for(i in 1:infos$n_z){
        Z_name = infos$n_z_vars[i]
        Z_pos = unique(stats::na.omit(infos$OUT$Pred_no[infos$OUT$Pred_Z == Z_name]))
        btw.Z[1:N,Z_pos] = stats::rnorm(n = N, mean = 0,
                                 sd = btw.var.sds[which(names(btw.var.sds) == Z_name)])
        #    colnames(btw.Z)[i] = Z_name
      }
    } else {
      btw.Z = btw_random
    }

    for(i in 1:infos$n_out){
      alpha = model$true.val[grepl(model$Param, pattern = paste0("alpha_",infos$out_var[i]))]
      sigma = model$true.val[grepl(model$Param, pattern = paste0("sigma_",infos$out_var[i]))]
      # create outcome values
      out_use = infos$OUT[infos$OUT$Var == infos$out_var[i],]
      if(nrow(out_use)>1){
        outs[,i] = alpha + btw.Z[,out_use$Pred_no]%*%out_use$true.val + stats::rnorm(n = N, mean = 0, sd = sigma)
      } else {
        outs[,i] = alpha + btw.Z[,out_use$Pred_no]*out_use$true.val + stats::rnorm(n = N, mean = 0, sd = sigma)
      }
    }

    OUT = data.frame(
      "ID" = 1:N,
      outs
    )
    colnames(OUT)[2:(infos$n_out+1)] = infos$out_var
  }



  # combine information =======================================================
  data = within

  if(infos$n_cov>1){
    W = cbind("ID" = 1:N, W)
    data = merge(x = data, y = W[,c("ID", infos$n_cov_vars)], by = "ID")
  }
  if(infos$n_out>0){
    data = merge(x = data, y = OUT, by = "ID")
    if(infos$n_z>0){
      colnames(btw.Z) = c(infos$re_pars$Param, infos$n_z_vars)
      btw.Z = cbind("ID" = 1:N, btw.Z)
      data = merge(x = data, y = btw.Z[,c("ID",infos$n_z_vars)])
    }
  }

  # return list
  VARsimData = list(
    model = model,
    data = data,
    RE.pars = btw_random
  )

  # add class
  class(VARsimData) <- "mlts_simdata"


  return(VARsimData)

}
