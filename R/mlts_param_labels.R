#' Convert parameter labels from mlts model to Stan labels and vice versa
#'
#' @param model `data frame`. Output of \code{\link[mlts]{mlts_model}} and
#' related functions.
#' @return An object of class `data.frame`.
#' @noRd
#'
#' @details
#' For internal use only.
#'
mlts_param_labels <- function(model){

  # eval model
  infos = mlts_model_eval(model)

  # add grouping variable if not present
  if(is.null(model$group)){
    model$group <- 1
  }

  # helper function to map names of parameters used in model with
  # parameter labels used in the stan models
  model$Param_stan = NA
#  model$Param_expr = model$Param
  model$Par_no = 1: nrow(model)

  # first add infos to model to extract variable names
  ##### FIXED EFFECT INTERCEPTS ================================================
  FEints = model[model$Type=="Fixed effect" & model$isRandom==1,]
  FEints$Param_stan = paste0("gammas[",FEints$group,",",rep(1:nrow(FEints[FEints$group==1,]),infos$G),"]")

  ##### CONSTANT DYNAMIC PARAMETERS ============================================
  FEdyn = model[model$Type=="Fixed effect" & model$isRandom==0 & model$Param_Label=="Dynamic",]
  if(nrow(FEdyn)>0){
    FEdyn$Param_stan = paste0("b_fix[",FEdyn$group,",",rep(1:nrow(FEdyn[FEdyn$group==1,]),infos$G),"]")
  }

  ##### CONSTANT INNOVATION VARIANCES ==========================================
  FEsigma = model[model$Type=="Fixed effect" & model$isRandom==0 & model$Param_Label=="Innovation Variance",]
  if(nrow(FEsigma)>0){
    FEsigma$Param_stan = paste0("sigma[",FEsigma$group,",",rep(1:nrow(FEsigma[FEsigma$group==1,]),infos$G),"]")
  }

  ##### RANDOM EFFECT SDs ======================================================
  REsds = model[model$Type == "Random effect SD",]
  REsds$Param_stan = paste0("sd_R[",REsds$group,",",rep(1:nrow(REsds[REsds$group==1,]),infos$G),"]")


  ##### RE CORRELATIONS ========================================================
  REcors = list()
  for(g in 1:infos$G){
    REcors[[g]] = model[model$Type == "RE correlation" & model$group==g,]
    if(nrow(REcors[[g]] > 0)){
      rand_pars = FEints$Param[FEints$group==g]
      rand_pars_pos = 1:length(rand_pars)
      REcors[[g]]$Param_stan = REcors[[g]]$Param
      for(i in 1:length(rand_pars)){
        REcors[[g]]$Param_stan = gsub(REcors[[g]]$Param_stan,
                               pattern = rand_pars[i],
                               replacement = rand_pars_pos[i], fixed = TRUE)
      }
      REcors[[g]]$Param_stan = gsub(REcors[[g]]$Param_stan, fixed = TRUE,
                               pattern = ".", replacement = ",")
      REcors[[g]]$Param_stan = gsub(REcors[[g]]$Param_stan, fixed = TRUE,
                                 pattern = "r_", replacement = paste0("bcorr[",g,","))
      REcors[[g]]$Param_stan = paste0(REcors[[g]]$Param_stan,"]")
    }
  }
  REcors = do.call(rbind,REcors)

  ###### INNOVATION COVARIANCE =================================================
  Fix.Covs <- list()
  for(g in 1:infos$G){
    Fix.Covs[[g]] = model[startsWith(model$Param, prefix = "r.zeta") & model$group==g,]

    if(nrow(Fix.Covs[[g]])>0){
      if(infos$D_cen == infos$q){
      Fix.Covs[[g]]$Param_stan = paste0("bcorr_inn[",g,",",
                                 substr(Fix.Covs[[g]]$Param, start = 8, 8),",",
                                 substr(Fix.Covs[[g]]$Param, start = 9, 9),"]")
      } else {
        D1 = as.integer(substr(Fix.Covs[[g]]$Param, start = 8, 8))
        D2 = as.integer(substr(Fix.Covs[[g]]$Param, start = 9, 9))
        # replace
        D1 = unlist(lapply(D1, function(x){infos$D_cen_pos[x]}))
        D2 = unlist(lapply(D2, function(x){infos$D_cen_pos[x]}))
        Fix.Covs[[g]]$Param_stan = paste0("bcorr_inn[",g,",",D1,",",D2,"]")
        }
    }
  }
  Fix.Covs = do.call(rbind, Fix.Covs)

  ###### RE on BETWEEN-LEVEL COVARIATES ========================================
  REpred <- list()
  for(g in 1:infos$G){
    REpred[[g]] = model[model$Type == "RE prediction" & model$group == g,]
    if(nrow(REpred[[g]])>0){
      infos$RE.PREDS$Param_stan = paste0("b_re_pred[",g,",",infos$RE.PREDS$re_pred_b_no,"]")
      REpred[[g]]$Param_stan <- NULL
      REpred[[g]] = merge(REpred[[g]],y = infos$RE.PREDS[,c("Param", "Param_stan")],
                   by = "Param", sort = FALSE)
      }
    }
  REpred <- do.call(rbind, REpred)


  ###### OUTCOME PREDICTION ====================================================
  # get the order of parameters in stan model

  OUTpred <- list()
  for(g in 1:infos$G){
  OUTpred[[g]] = model[model$Type == "Outcome prediction" & model$group==g,]
  if(nrow(OUTpred[[g]]) > 0){
    # regression parameter
    infos$OUT = infos$OUT[order(infos$OUT$out_var_no, infos$OUT$Pred_no),]
    infos$OUT$Param_stan = paste0("b_out_pred[",g,",",1:nrow(infos$OUT),"]")

    # add the helper columns
    OUTpred[[g]]$Param_stan <- NULL
    OUTpred[[g]] = merge(OUTpred[[g]], y = infos$OUT[,c("Param", "out_var_no", "Param_stan")],
                    by = c("Param"), all = T, sort = FALSE)

    for(i in 1:infos$n_out){
      OUTpred[[g]]$out_var_no[endsWith(OUTpred[[g]]$Param, infos$out_var[i])] = i
    }
      # intercepts and residual SDs
    OUTpred[[g]]$Param_stan[OUTpred[[g]]$Param_Label == "intercept"] <- paste0(
      "alpha_out[",g,",",OUTpred[[g]]$out_var_no[OUTpred[[g]]$Param_Label == "intercept"],"]")
    OUTpred[[g]]$Param_stan[OUTpred[[g]]$Param_Label == "Residual SD"] <- paste0(
      "sigma_out[",g,",",OUTpred[[g]]$out_var_no[OUTpred[[g]]$Param_Label == "Residual SD"],"]")

    # remove helper columns
    OUTpred[[g]] = OUTpred[[g]][,colnames(FEints)]
    }
  }
  OUTpred = do.call(rbind, OUTpred)


  ### MEASUREMENT MODEL PARAMETERS =============================================
  if(infos$isLatent == TRUE){
    N_inds = nrow(infos$indicators)
    alphas = data.frame(
      "Param" = paste0("alpha_",infos$indicators$q, ".",infos$indicators$p),
      "Param_stan" = paste0("alpha[",1:N_inds,"]")#,
      )
    loadB = data.frame(
      "Param" = paste0("lambdaB_",infos$indicators$q, ".",infos$indicators$p),
      "Param_stan" = paste0("loadB[",1:N_inds,"]")#,
      )
    loadW = data.frame(
      "Param" = paste0("lambdaW_",infos$indicators$q, ".",infos$indicators$p),
      "Param_stan" = paste0("loadW[",1:N_inds,"]")#,
      )
    sigmaB = data.frame(
      "Param" = paste0("sigmaB_",infos$indicators$q, ".",infos$indicators$p),
      "Param_stan" = paste0("sigmaB[",1:N_inds,"]")#,
      )
    sigmaW = data.frame(
      "Param" = paste0("sigmaW_",infos$indicators$q, ".",infos$indicators$p),
      "Param_stan" = paste0("sigmaW[",1:N_inds,"]")#,
      )
    mm.pars = rbind(alphas, loadB, sigmaB, loadW, sigmaW)


    mm.pars = merge(x = mm.pars,
                    y = model[model$Model=="Measurement",colnames(model) != "Param_stan"],
                    by = "Param", all.y = T)

  }

  #### COMBINE
  par_tab = rbind(FEints, FEdyn, Fix.Covs, FEsigma, REsds, REcors, REpred, OUTpred)

  if(infos$isLatent == TRUE){
    # par_tab = plyr::rbind.fill(par_tab, mm.pars)
    par_tab <- dplyr::bind_rows(par_tab, mm.pars)
  }

  return(par_tab)
}



