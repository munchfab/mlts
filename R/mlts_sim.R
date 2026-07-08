#' Simulate data from mlts model
#'
#' @details
#' A function to generate data from an output of \code{\link[mlts]{mlts_model}}.
#'
#' @param model `data.frame`. Output of \code{\link[mlts]{mlts_model}}.
#' @param default logical. If set to `TRUE`, default prior specifications are
#' added.
#' @param N integer Number of observational units.
#' @param N_G vector of integers. Number of observational units per group.
#' @param TP integer. Number of measurements per observational unit.
#' @param burn.in integer. Length of ‘burn-in’ period.
#' @param seed integer. Seed used for data generation.
#' @param seed.true integer. Separate seed used for sampling of true
#' population parameters values from plausible ranges for stationary time series.
#' @param btw.var.sds named numeric vector. Provide standard deviation(s) for all exogenous
#' between-level variable(s) specified in `model`, e.g. (`btw.var.sds = c("covariate1" = 1)`,
#' to set the SD of the variable "covariate1" to 1). Mean values of the respective
#' variable(s) will be set to 0 per default.
#' @param exogenous Matrix of numeric values of exogenous variables with `N`*(`TP`+`burn.in`)
#' rows and separate columns for each variable.
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

mlts_sim <- function(model, default = FALSE, N = NULL, N_G = NULL, TP, burn.in = 50, seed = NULL,
                     seed.true = 1, btw.var.sds = NULL, exogenous = NULL){



  # use helper function to read out information on model
  infos = mlts_model_eval(model)

  # check correct inputs:
  if(infos$G > 1 & is.null(N_G)){
    stop("Use `N_G` to specify the number of clusters in each group.")
  }
  if(infos$G > 1 & !is.null(N_G)){
    N <- sum(N_G)
  }

  # set seed
  set.seed(seed.true)

  # check that exogenous variable is entered and in the correct format
  if(any(infos$is_wcen==0)){

    if(is.null(exogenous)){
      stop("Enter values of exogenous variables as a matrix with the correct dimensions:
          rows = N*(TP+burn.in)
          columns = no of exogenous variables")
    }
    check2 = dim(exogenous)[1] == N*(TP+burn.in)
    check3 = dim(exogenous)[2] == sum(infos$is_wcen == 0)

    if(check2 == FALSE | check3 == FALSE){
      stop("Enter values of exogenous variables as a matrix with the correct dimensions:
          rows = N*(TP+burn.in)
          columns = no of exogenous variables")
    }
  }


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
    model$true.val[model$Level == "Within" & model$Type == "Loading" & model$Constraint != "= 1"] =
      sample(x = seq(0.7, 0.9, by = 0.05), size = infos$n_loadWfree, replace = TRUE)
    model$true.val[model$Level == "Between" & model$Type == "Loading" & model$Constraint == "= 1"] = 1
    model$true.val[model$Level == "Between" & model$Type == "Loading" & model$Constraint != "= 1"] =
      sample(x = seq(0.7, 0.9, by = 0.05), size = infos$n_loadB, replace = TRUE)

    # evaluate equality constraints on loading parameters across levels
    equal = unique(model$Constraint[model$Type == "Loading" & !(model$Constraint %in% c("= 1", "free"))])
    n_equal = length(equal)
    if(n_equal>0){
      for(z in equal){
        # set all loadings with the same equality constraint to the same value
        val = model$true.val[model$Constraint %in% z]
        model$true.val[model$Constraint %in% z] <- val[1]
      }
    }

    model$true.val[model$Level == "Within" & model$Type == "Measurement Error SD" & model$Constraint == "= 0"] = 0
    model$true.val[model$Level == "Within" & model$Type == "Measurement Error SD" & model$Constraint == "free"] =
      sample(x = seq(0.15, 0.3, by = 0.05), size = infos$n_sigmaWfree, replace = TRUE)
    model$true.val[model$Level == "Between" & model$Type == "Measurement Error SD" & model$Constraint == "= 0"] = 0
    model$true.val[model$Level == "Between" & model$Type == "Measurement Error SD" & model$Constraint == "free"] =
      sample(x = seq(0.15, 0.2, by = 0.05), size = infos$n_sigmaBfree, replace = TRUE)


    # DYNAMIC PART =============================================================
    ## add helper column, in case no group is specified
    if( infos$G == 1 ){ model$group <- 1 }

    for ( gg in 1:infos$G ){

      # FIXED EFFECTS:
      type = "Fixed effect"

      ## MUS =====
      use      <- startsWith(model$Param_Label, "Trait")
      vals     <- sample(x = seq(from= 0, to = 1, by = 0.1), size = 1000, replace = T)
      model    <- add_trues(model, type = type, group = gg, which = use, values = vals, adjust_size = TRUE)
      # --

      ## PHIs ====

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
            phis$true.val[i] = phi.lag1/as.numeric(phis$Lag[i])
          }
        }
      }

      vals  <- round(phis$true.val,3)
      model <- add_trues(model, type = type, label = "Dynamic", group = gg, values = vals)
      #-

      ## t0-effect
      t0_effect = infos$fix_pars_dyn[startsWith(infos$fix_pars_dyn$Param,"phi(0"),]
      if(nrow(t0_effect) > 0){
        for(i in 1:nrow(t0_effect)){
          # choose value according to the respective lagged effect - if present
          lag_value = model$true.val[model$group == gg, model$Param == paste0("phi(1)_",t0_effect$Dout[i],t0_effect$Dpred[i])]

          if(is.numeric(lag_value)){
            vals <- 1.2 * lag_value
          } else {
            vals <- sample(x = seq(from=-0.15,to=0.15,by=0.05),size = 1)
          }

          use = model$Param == t0_effect$Param[i]
          model <- add_trues(model, group = gg, which = use, values = vals)
        }
      }
      # -

      ## interaction effects
      if(infos$n_int > 0){
        use  <- startsWith(model$Param, prefix = "phi(i)")
        vals <- sample(x = seq(from=-0.1,to=0.1,by=0.025), size = infos$n_int)
        model <- add_trues(model, group = gg, which = use, values = vals)
      }
      # --

      ## LOG INNOVATION (CO)VARIANCES
      model <- add_trues(model, type = type, label = "Log Innovation Variance", group = gg, values = -0.3)
      model <- add_trues(model, type = type, label = "Innovation Variance", group = gg, values = 0.75)
      model <- add_trues(model, type = type, label = "Log Innovation Covariance", group = gg, values = 0.3)
      model <- add_trues(model, type = type, label = "Log Innovation Covariance", group = gg, values = -0.3)
      model <- add_trues(model, type = type, label = "Innovation correlation", group = gg, values = -0.15)

      # ---


      # RANDOM EFFECTS SDs:
      type = "Random effect SD"

      ## Mus
      use   <- startsWith(model$Param_Label, "Trait")
      vals  <- sample(seq(from= 0.7, to = 1.2, by = 0.1), size = 100, replace = TRUE)
      model <- add_trues(model, type = type, which = use, group = gg, values = vals, adjust_size = TRUE)

      ## Phis
      model <- add_trues(model, type = type, label = "Dynamic", group = gg, values = 0.15)
      ## smaller for interaction effects
      use = startsWith(model$Param, "phi(i)")
      model <- add_trues(model, type = type, which = use, group = gg, values = 0.1)
      ## log innovation variances
      model <- add_trues(model, type = type, label = "Log Innovation Variance", group = gg, values = 0.25)
      ## log innovation covariance(s)
      model <- add_trues(model, type = type, label = "Log Innovation Covariance", group = gg, values = 0.25)

      # ---

      # RANDOM EFFECTS CORRELATIONS:

      ## set all to zero for now
      type <- "RE correlation"
      vals <- sample(seq(from = -0.2, to = 0.2, by = 0.05), replace = TRUE, size = sum(model$Type == type & model$group==gg))
      model <- add_trues(model = model, type = type, group = gg, values = vals)



      # ---

      # RE as OUTCOME:
      type <- "RE prediction"
      vals <- sample(seq(from = -0.2, to = 0.2, by = 0.05), replace = TRUE, size = sum(model$Type == type & model$group==gg))
      model <- add_trues(model, type = type, group = gg, values = vals)

      # ---

      # OUTCOME PREDICTION :
      type <- "Outcome prediction"
      vals <- sample(seq(from = -0.3, to = 0.3, by = 0.1), replace = TRUE, size = sum(model$Type == type & model$group==gg))
      model <- add_trues(model, type = type, group = gg, values = vals)

      # scale true values for AR and CL as predictor
      model$true.val[model$Type == type & grepl(pattern = "phi",model$Param)] <-
        model$true.val[model$Type == type & grepl(pattern = "phi",model$Param)] * 5

      model <- add_trues(model, type = type, label ="intercept", group = gg, values = 0)
      model <- add_trues(model, type = type, label ="Residual SD", group = gg, values = 0.5)

      # ---

      }

      model$true.val <- round(model$true.val, 3)

      } else if ( is.null(model$true.val)) {
      stop(
        paste("No true parameter values provided in model$true.val. Set default = TRUE to run data generation with random true parameter values.",
         "Alternatively, user-specified values for each parameter can be specified in an additional column `true.val` in model.")
      )
    }



  # set seed for data generation
  if(!is.null(seed)){
    set.seed(seed)
  }

  # run again after adding true parameter values
  if(infos$G == 1){
    infos = mlts_model_eval(model)
  } else {
    infos = mlts_model_eval(model)
  }



  # start generating between-level model =======================================

  # separately by group
  if( infos$G == 1 ) { N_G <- N }

  # store final person parameters as a list
  btw <- list()
  W <- list()
  for( gg in 1:infos$G ){

  # FIXED EFFECTS
  gammas = model$true.val[model$group == gg & model$Type=="Fixed effect" & model$isRandom==1]

  # BETWEEN-LEVEL
  # sample covariates and get expected values of individual parameters
  bmu = matrix(data = NA, nrow = N_G[gg], ncol = infos$n_random)
  W[[gg]] = matrix(data = NA, nrow = N_G[gg], ncol = infos$n_cov)
  cov_name = c()
  W[[gg]][,1] = 1 # intercept
  if(infos$n_cov>1){
    for(i in 2:infos$n_cov){
      cov_name[i-1] = unique(infos$RE.PREDS$re_preds[infos$RE.PREDS$pred_no == i-1])
      W[[gg]][,i] = stats::rnorm(n = N_G[gg], mean = 0, btw.var.sds[names(btw.var.sds) == cov_name[i-1]])
    }
  }
  colnames(W[[gg]]) <- c("Intercept", cov_name)

  for(i in 1:infos$n_random){
    # get expected individual parameters
    pred_use = infos$RE.PREDS[infos$RE.PREDS$re_no ==i,]
    val_use = model$true.val[model$group==gg & model$Param %in% pred_use$Param]

    if(nrow(pred_use)>0){
      bmu[,i] = W[[gg]][,c(1,pred_use$pred_no+1)] %*% c(gammas[i], val_use)
    } else {
      bmu[,i] = gammas[i]
    }
  }

  # calculate covariances from correlations
  n_random = infos$n_random
  rand.pars = infos$re_pars$Param

  # variance covariance matrix of random effects
  if(n_random == 1){
    cov_mat = model$true.val[model$Type=="Random effect SD" & model$group == gg]
  } else {
    cov_mat = diag(model$true.val[model$Type=="Random effect SD" & model$group == gg]^2)
    for(i in 1:n_random){
      for(j in 1:n_random){
        if(i < j){
          r = model$true.val[model$Param == paste0("r_",rand.pars[i],".", rand.pars[j]) & model$group == gg]
          cov_mat[i,j] = cov_mat[j,i] <- r * sqrt(cov_mat[i,i]) * sqrt(cov_mat[j,j])
        }
      }
    }
  }


  #### sample random effects from multivariate normal distribution and add to bmus
  btw_random = matrix(NA, nrow = N_G[gg], ncol = infos$n_random)
  if(n_random == 1){
    btw_random = bmu + stats::rnorm(n = N_G[gg], mean = 0, sd = cov_mat)
  } else {
    btw_random = bmu + mvtnorm::rmvnorm(n = N_G[gg], mean = rep(0, infos$n_random), sigma = cov_mat)
  }
  colnames(btw_random) = infos$fix_pars$Param[infos$is_random]

  # check for AR parameters with absolute values below "1"
  posAR = infos$fix_pars[infos$fix_pars$isRandom==1 & infos$fix_pars$isAR==1, "no"]

  if(sum(abs(btw_random[,posAR]) >= 1)){
    stop("Absolute individual AR effects greater than 1 were sampled. Consider
            setting the true values of the fixed effect or the random effect SD to a lower value.")
  }

  # now combine fixed effects and random effects
  btw[[gg]] = matrix(NA, nrow = N_G[gg], infos$n_pars)
  btw[[gg]][, infos$is_random] = btw_random               # first add random pars
  if(infos$n_fixed>0){
    btw[[gg]][,infos$is_fixed[1,]] = rep(model$true.val[model$Type=="Fixed effect" & model$group == gg][infos$is_fixed[1,]], each=N_G[gg])
  }
  if(infos$n_innos_fix>0){
    for(i in infos$innos_fix_pos)
      btw[[gg]][,infos$innos_pos[i]] = rep(model$true.val[model$Type=="Fixed effect" & model$Param_Label == "Innovation Variance"& model$group == gg][i],times=N_G[gg])
    }
  }

  #### WITHIN-LEVEL PROCESS ====================================================
  btw <- do.call(rbind, btw)
  mm_pars = model
  cor_pars = model
  mm_pars$sample <- mm_pars$true.val
  cor_pars$sample <- cor_pars$true.val

  # create vector of group_ids
  group_ids <- unlist(sapply(1:infos$G, function(x) rep(x, N_G[x])))

  within <- mlts_sim_within(
    infos = infos,
    N = sum(N_G),
    group_ids = group_ids,
    TP = TP,
    burn.in = burn.in,
    btw = btw,
    mm_pars = mm_pars,
    cor_pars = cor_pars,
    exogenous = exogenous
  )

    # CREATE OUTCOMES ==========================================================
    outs = matrix(NA, ncol = infos$n_out, nrow = N)
    if(infos$n_out > 0){
      if(infos$n_z > 0){
        btw.Z = matrix(nrow = N, ncol = (infos$n_random + infos$n_z))
        btw.Z[,1:infos$n_random] = btw
        for(i in 1:infos$n_z){
          Z_name = infos$n_z_vars[i]
          Z_pos = unique(stats::na.omit(infos$OUT$Pred_no[infos$OUT$Pred_Z == Z_name]))
          btw.Z[1:N,Z_pos] = stats::rnorm(n = sum(N_G), mean = 0,
                                          sd = btw.var.sds[which(names(btw.var.sds) == Z_name)])
          #    colnames(btw.Z)[i] = Z_name
        }
      } else {
        btw.Z = btw
      }

      for(i in 1:infos$n_out){
        for(gg in 1:infos$G){
          alpha = model$true.val[model$group == gg & grepl(model$Param, pattern = paste0("alpha_",infos$out_var[i]))]
          sigma = model$true.val[model$group == gg & grepl(model$Param, pattern = paste0("sigma_",infos$out_var[i]))]
          # create outcome values
          out_use = infos$OUT[infos$OUT$Var == infos$out_var[i],]
          if(nrow(out_use)>1){
            outs[group_ids==gg,i] = alpha + btw.Z[group_ids==gg,out_use$Pred_no]%*%out_use$true.val + stats::rnorm(n = N_G[gg], mean = 0, sd = sigma)
          } else {
            outs[group_ids==gg,i] = alpha + btw.Z[group_ids==gg,out_use$Pred_no]*out_use$true.val + stats::rnorm(n = N_G[gg], mean = 0, sd = sigma)
          }
        }
      }

      OUT = data.frame(
        "ID" = 1:N,
        outs
      )
      colnames(OUT)[2:(infos$n_out+1)] = infos$out_var
    }


    # add censored versions of variables if requested
    y_cols = colnames(within)[startsWith(colnames(within),"Y")]
    if(!is.null(attr(model, which = "censor_left"))){

      # get censoring threshold
      censor_LB <- attr(model, which = "censor_left")

      within[,paste0(y_cols,"_cens")] <- within[,y_cols]
      within[,paste0(y_cols,"_cens")][within[,y_cols] <= censor_LB] <- censor_LB
    }

    if(!is.null(attr(model, which = "censor_right"))){

      # get censoring threshold
      censor_UB <- attr(model, which = "censor_right")
      if(is.null(attr(model, which = "censor_left"))){
        within[,paste0(y_cols,"_cens")] <- within[,y_cols]
      }
      within[,paste0(y_cols,"_cens")][within[,paste0(y_cols,"_cens")] >= censor_UB] <- censor_UB
    }

    # combine information =======================================================
    data = within

    if(infos$n_cov>1){
      W = cbind("ID" = 1:N, do.call(rbind, W))
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
