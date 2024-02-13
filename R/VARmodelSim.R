#' Title
#'
#' @param VARmodel data.frame. Output of VARmodel-Functions.
#' @param default logical. If set to `TRUE`, default prior specifications are
#' added.
#'
#' @return An object of class `data.frame`.
#' @export
#'
VARmodelSim <- function(VARmodel, default = F, N, TP, burn.in = 500, seed = NULL,
                        btw.var.sds = NULL){

  if(!is.null(seed)){
    set.seed(seed)
  }

  ### start with a check that sds are entered for all between variables


  # simulate data based on VARmodel

  ### write separate function for specification of true parameter values
  ### for now use fixed values


  # use helper function to read out information on model
  infos = VARmodelEval(VARmodel)


  # use some default settings for parameter values
  if(default==T){
    VARmodel$true.val = NA

    # MEASUREMENT MODEL PARAMETERS =============================================
    n_p = sum(infos$p)
    model = "Measurement"
    VARmodel$true.val[VARmodel$Model == model & VARmodel$Type == "Item intercepts" & VARmodel$Constraint == "= 0"] = 0
    VARmodel$true.val[VARmodel$Model == model & VARmodel$Type == "Item intercepts" & VARmodel$Constraint == "free"] =
      sample(x = seq(0.5, 2, by = 0.5), size = infos$n_alphafree, replace = T)
    VARmodel$true.val[VARmodel$Level == "Within" & VARmodel$Type == "Loading" & VARmodel$Constraint == "= 1"] = 1
    VARmodel$true.val[VARmodel$Level == "Within" & VARmodel$Type == "Loading" & VARmodel$Constraint == "free"] =
      sample(x = seq(0.7, 0.9, by = 0.05), size = infos$n_loadWfree, replace = T)
    VARmodel$true.val[VARmodel$Level == "Between" & VARmodel$Type == "Loading" & VARmodel$Constraint == "= 1"] = 1
    VARmodel$true.val[VARmodel$Level == "Between" & VARmodel$Type == "Loading" & VARmodel$Constraint == "free"] =
      sample(x = seq(0.7, 0.9, by = 0.05), size = infos$n_loadBfree, replace = T)
    VARmodel$true.val[VARmodel$Level == "Within" & VARmodel$Type == "Measurement Error SD" & VARmodel$Constraint == "= 0"] = 0
    VARmodel$true.val[VARmodel$Level == "Within" & VARmodel$Type == "Measurement Error SD" & VARmodel$Constraint == "free"] =
      sample(x = seq(0.2, 0.3, by = 0.05), size = infos$n_sigmaWfree, replace = T)
    VARmodel$true.val[VARmodel$Level == "Between" & VARmodel$Type == "Measurement Error SD" & VARmodel$Constraint == "= 0"] = 0
    VARmodel$true.val[VARmodel$Level == "Between" & VARmodel$Type == "Measurement Error SD" & VARmodel$Constraint == "free"] =
      sample(x = seq(0.2, 0.3, by = 0.05), size = infos$n_sigmaBfree, replace = T)

    # FIX EFFECTS ========
    model.type = "Fix effect"
    ## Mus
    VARmodel$true.val[VARmodel$Type==model.type & startsWith(VARmodel$Param_Label, "Trait")] = 1
    ## Phis
    VARmodel$true.val[VARmodel$Type==model.type & VARmodel$Param_Label=="Dynamic"] = sample(
      c(round(seq(from = 0.1, to = 0.3, by = 0.05),3)), replace = T,
      size = nrow(infos$fix_pars_dyn))
    ## log innovation variances
    VARmodel$true.val[VARmodel$Type==model.type & VARmodel$Param_Label=="Log Innovation Variance"] = -0.3
    ## Fixed innovation variance
    VARmodel$true.val[VARmodel$Type==model.type & VARmodel$Param_Label=="Innovation Variance"] = 0.75

    ###### NEEDS UPDATING ----
    ## Log innovation covariance
    VARmodel$true.val[VARmodel$Type==model.type & VARmodel$Param_Label=="Log Innovation Covariance"] = -1
    ###### ----


    # RANDOM EFFECT SDs ==========
    model.type = "Random effect SD"
    ## Mus
    VARmodel$true.val[VARmodel$Type==model.type & startsWith(VARmodel$Param_Label, "Trait")] = 1
    ## Phis
    VARmodel$true.val[VARmodel$Type==model.type & VARmodel$Param_Label=="Dynamic"] = 0.15
    ## log innovation variances
    VARmodel$true.val[VARmodel$Type==model.type & VARmodel$Param_Label=="Log Innovation Variance"] = 0.4
    ## log innovation covaraince(s)
    VARmodel$true.val[VARmodel$Type==model.type & VARmodel$Param_Label=="Log Innovation Covariance"] = 0.15


    # RANDOM EFFECT CORRELATIONS ============
    ## set all to zero for now
    model.type = "RE correlation"
    VARmodel$true.val[VARmodel$Type == model.type] = sample(
      c(round(seq(from = -0.3, to = 0.3, by = 0.1),3)), replace = T,
      size = sum(VARmodel$Type == model.type))

    # RE as OUTCOME =========================
    model.type = "RE prediction"
    VARmodel$true.val[VARmodel$Type == model.type] = sample(
      c(round(seq(from = -0.2, to = 0.2, by = 0.05),3)), replace = T,
        size = sum(VARmodel$Type == model.type))

    # OUTCOME PREDICTION ===================
    model.type = "Outcome prediction"
    VARmodel$true.val[VARmodel$Type == model.type] = sample(
      c(round(seq(from = -0.3, to = 0.3, by = 0.1),3)), replace = T,
      size = sum(VARmodel$Type == model.type))
    VARmodel$true.val[VARmodel$Type == model.type & VARmodel$Param_Label == "intercept"] = 0
    VARmodel$true.val[VARmodel$Type == model.type & VARmodel$Param_Label == "Residual SD"] = 0.5

  }

  # run again after adding true parameter values
  infos = VARmodelEval(VARmodel)


  # start generating between-level model =======================================

  # FIXED EFFECTS
  gammas = VARmodel$true.val[VARmodel$Type=="Fix effect" & VARmodel$isRandom==1]

  # BETWEEN-LEVEL
  # sample covariates and get expected values of individual parameters
  bmu = matrix(data = NA, nrow = N, ncol = infos$n_random)
  W = matrix(data = NA, nrow = N, ncol = infos$n_cov)
  cov_name = c()
  W[,1] = 1 # intercept
  if(infos$n_cov>1){
    for(i in 2:infos$n_cov){
      cov_name[i-1] = unique(infos$RE.PREDS$re_preds[infos$RE.PREDS$pred_no == i-1])
      W[1:N,i] = rnorm(n = N, mean = 0, btw.var.sds[names(btw.var.sds) == cov_name[i-1]])
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

  # variance covariance matrix of random effects
  cov_mat = diag(VARmodel$true.val[VARmodel$Type=="Random effect SD"]^2)

  # calculate covariances from correlations
  n_random = infos$n_random
  rand.pars = infos$re_pars$Param
  for(i in 1:n_random){
    for(j in 1:n_random){
      if(i < j){
        r = VARmodel$true.val[VARmodel$Param == paste0("r_",rand.pars[i],".", rand.pars[j])]
        cov_mat[i,j] = cov_mat[j,i] <- r * sqrt(cov_mat[i,i]) * sqrt(cov_mat[j,j])
      }
    }
  }


  #### sample random effects from multivariate normal distribution and add to bmus
  btw_random = matrix(NA, nrow = N, ncol = infos$n_random)
  btw_random = bmu + mvtnorm::rmvnorm(n = N, mean = rep(0, infos$n_random), sigma = cov_mat)
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
    btw[,infos$is_fixed] = VARmodel$true.val[VARmodel$Type=="Fix effect"][infos$is_fixed[1,]]
  }
  if(infos$n_innos_fix>0){
    for(i in infos$innos_fix_pos)
    btw[,infos$innos_pos[i]] = VARmodel$true.val[VARmodel$Type=="Fix effect" & VARmodel$Param_Label == "Innovation Variance"][i]
  }

  #### WITHIN-LEVEL PROCESS ====================================================
  # prepare data frame
  NT = TP + burn.in
  within = data.frame(
    "ID" = rep(1:N, each = NT),
    "time" = rep(1:NT, times = N)
  )
  q = infos$n_mus   # number of constructs
  y_cols = paste0("Y",1:q)   # prepare columns
  within[, y_cols] = NA

  for(i in 1:N){
    # build person-specific transition matrix
    transition = matrix(nrow = q, ncol = q, data = 0)
    # build person-specific prediction error matrix
    innoVars.i = btw[i,infos$innos_pos]
    for(d in 1:q){ # loop over number of constructs
      transition[d,] = btw[i,infos$Dpos1[d]:infos$Dpos2[d]]
      if(infos$innos_rand[d] == 1){
        innoVars.i[d] = exp(innoVars.i[d]) # retransform log innovations
      } else{
        innoVars.i[d] = innoVars.i[d]^2
      }
    }

    if(q == 1){
      inno_var_mat = matrix(data = innoVars.i, nrow = 1, ncol = 1)
    } else {
      inno_var_mat = diag(innoVars.i)

      ########## ADD INNOVATION COVARIANCES HERE ::::::::::::::::::::::::::::::::



      # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    }

    for(t in 1:NT){
      # use means as starting values
      if(within[within$ID==i & within$time == t,"time"] == 1){
        within[within$ID==i & within$time == t,y_cols] = btw[i,1:q]
      } else {
        y_lag = within[within$ID==i & within$time == t-1,y_cols]
        if(q > 1){
         y = as.matrix(y_lag) %*% t(transition) # get expected values
        } else {
         y = y_lag * transition
        }
        y = y + mvtnorm::rmvnorm(n = 1, mean = rep(0, q), sigma = inno_var_mat) # add error
        within[within$ID==i & within$time==t, y_cols] = y
      }
    }

    # add trait scores (for manifest indicators)
    if(infos$isLatent == F){
      for(j in 1:infos$n_mus){
      within[within$ID==i,y_cols[j]] = btw[i,j] + within[within$ID==i,y_cols[j]]
      }
    }
  }

  # remove burn.in
  within$time = within$time - burn.in
  within = within[within$time>0,]
  # --------


  ##### add measurement model here -------------------------------------------
  # create manifest indicator scores
  if(infos$isLatent == TRUE){
  N_inds = max(infos$indicators$p_pos)
  for(i in 1:N_inds){
    q = infos$indicators$q[i]
    p = infos$indicators$p[i]
    ind.lab = paste0("Y",q,".",p)
    # within
    loadW = VARmodel$true.val[VARmodel$Level == "Within" & VARmodel$Type == "Loading"][i]
    sigmaW = VARmodel$true.val[VARmodel$Level == "Within" & VARmodel$Type == "Measurement Error SD"][i]
    # between
    alpha = VARmodel$true.val[VARmodel$Param == paste0("alpha_",q,".",p)]
    alpha = ifelse(length(alpha) == 0, 0, alpha)
    loadB = VARmodel$true.val[VARmodel$Param == paste0("lambdaB_",q,".",p)]
    loadB = ifelse(length(loadB) == 0, 1, loadB)
    sigmaB = VARmodel$true.val[VARmodel$Param == paste0("sigmaB_",q,".",p)]
    sigmaB = ifelse(length(sigmaB) == 0, 0, sigmaB)

    for(p in 1:N){
      # create WITHIN-PART:
      etaW = within[within$ID==p, paste0("Y",q)]
      YW = loadW * etaW + ifelse(sigmaW == 0, 0, rnorm(n = N, mean = 0, sd = sigmaW))

      # create BETWEEN-PART:
      YB = alpha + loadB * btw[p,infos$indicators$etaB_pos[i]] + ifelse(sigmaB == 0, 0, rnorm(n = 1, mean = 0, sd = sigmaB))

      # indicator
      within[within$ID==p,ind.lab] = YB + YW
    }
  }
  }

  # CREATE OUTCOMES ==========================================================
  outs = matrix(NA, ncol = infos$n_out, nrow = N)
  if(infos$n_out > 0){
    if(infos$n_z > 0){
      btw.Z = matrix(nrow = N, ncol = (infos$n_random + infos$n_z))
      btw.Z[,1:infos$n_random] = btw_random
      for(i in 1:infos$n_z){
        Z_name = infos$n_z_vars[i]
        Z_pos = unique(na.omit(infos$OUT$Pred_no[infos$OUT$Pred_Z == Z_name]))
        btw.Z[1:N,Z_pos] = rnorm(n = N, mean = 0,
                                 sd = btw.var.sds[which(names(btw.var.sds) == Z_name)])
    #    colnames(btw.Z)[i] = Z_name
      }
    } else {
      btw.Z = btw_random
    }

    for(i in 1:infos$n_out){
      alpha = VARmodel$true.val[grepl(VARmodel$Param, pattern = paste0("alpha_",infos$out_var[i]))]
      sigma = VARmodel$true.val[grepl(VARmodel$Param, pattern = paste0("sigma_",infos$out_var[i]))]
      # create outcome values
      out_use = infos$OUT[infos$OUT$Var == infos$out_var[i],]
      if(nrow(out_use)>1){
         outs[,i] = alpha + btw.Z[,out_use$Pred_no]%*%out_use$true.val + rnorm(n = N, mean = 0, sd = sigma)
       } else {
         outs[,i] = alpha + btw.Z[,out_use$Pred_no]*out_use$true.val + rnorm(n = N, mean = 0, sd = sigma)
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
   VARmodel = VARmodel,
   data = data,
   RE.pars = btw_random
 )

 # add class
 class(VARsimData) <- "VARsimData"


 return(VARsimData)

}
