#' Title
#'
#' @param VARmodel data frame.
#' @export
#'
VARmodelParLabels <- function(VARmodel){

  # eval model
  infos = VARmodelEval(VARmodel)

  # helper function to map names of parameters used in VARmodel with
  # parameter labels used in the stan models
  VARmodel$Param_stan = NA
  VARmodel$Par_no = 1: nrow(VARmodel)

  # first add infos to VARmodel to extract variable names
  ##### FIXED EFFECT INTERCEPTS
  FEints = VARmodel[VARmodel$Type=="Fix effect" & VARmodel$isRandom==1,]
  FEints$Param_stan = paste0("gammas[",1:nrow(FEints),"]")

  ##### CONSTANT DYNAMIC PARAMETERS
  FEdyn = VARmodel[VARmodel$Type=="Fix effect" & VARmodel$isRandom==0 & VARmodel$Param_Label=="Dynamic",]
  if(nrow(FEdyn)>0){
  FEdyn$Param_stan = paste0("b_fix[",1:nrow(FEdyn),"]")
  }

  ##### CONSTANT INNOVATION VARIANCES
  FEsigma = VARmodel[VARmodel$Type=="Fix effect" & VARmodel$isRandom==0 & VARmodel$Param_Label=="Innovation Variance",]
  if(nrow(FEsigma)>0){
    FEsigma$Param_stan = paste0("sigma[",1:nrow(FEsigma),"]")
  }

  ##### RANDOM EFFECT SDs
  REsds = VARmodel[VARmodel$Type == "Random effect SD",]
  REsds$Param_stan = paste0("sd_R[",1:nrow(REsds),"]")


  ##### RE CORRELATIONS
  REcors = VARmodel[VARmodel$Type == "RE correlation",]
  if(nrow(REcors > 0)){
    rand_pars = FEints$Param
    rand_pars_pos = 1:length(rand_pars)
    REcors$Param_stan = REcors$Param
    for(i in 1:length(rand_pars)){
      REcors$Param_stan = gsub(REcors$Param_stan,
                                   pattern = rand_pars[i],
                                   replacement = rand_pars_pos[i])
    }
    REcors$Param_stan = gsub(REcors$Param_stan, fixed = T,
                                 patter = ".", replacement = ",")
    REcors$Param_stan = gsub(REcors$Param_stan, fixed = T,
                                 patter = "r_", replacement = "bcorr[")
    REcors$Param_stan = paste0(REcors$Param_stan,"]")
  }


  ###### RE on BETWEEN-LEVEL COVARIATES
  REpred = VARmodel[VARmodel$Type == "RE prediction",]
  if(nrow(REpred)>0){
    infos$RE.PREDS$Param_stan = paste0("b_re_pred[",infos$RE.PREDS$re_pred_b_no,"]")
    REpred$Param_stan <- NULL
    REpred = merge(REpred,y = infos$RE.PREDS[,c("Param", "Param_stan")],
                   by = "Param", sort = F)
  }

  ###### OUTCOME PREDICTION
  # get the order of parameters in stan model
  OUTpred = VARmodel[VARmodel$Type == "Outcome prediction",]
  if(nrow(OUTpred) > 0){
    # regression parameter
    infos$OUT = infos$OUT[order(infos$OUT$out_var_no, infos$OUT$Pred_no),]
    infos$OUT$Param_stan = paste0("b_out_pred[",1:nrow(infos$OUT),"]")

    # add the helper columns
    OUTpred$Param_stan <- NULL
    OUTpred = merge(OUTpred, y = infos$OUT[,c("Param", "out_var_no", "Param_stan")],
                    by = c("Param"), all = T, sort = F)

    for(i in 1:infos$n_out){
      OUTpred$out_var_no[endsWith(OUTpred$Param, infos$out_var[i])] = i
    }
      # intercepts and residual SDs
    OUTpred$Param_stan[OUTpred$Param_Label == "intercept"] <- paste0(
      "alpha_out[",OUTpred$out_var_no[OUTpred$Param_Label == "intercept"],"]")
    OUTpred$Param_stan[OUTpred$Param_Label == "Residual SD"] <- paste0(
      "sigma_out[",OUTpred$out_var_no[OUTpred$Param_Label == "Residual SD"],"]")

    # remove helper columns
    OUTpred = OUTpred[,colnames(FEints)]
  }

  #### COMBINE
  par_tab = rbind(FEints, FEdyn, FEsigma, REsds, REcors, REpred, OUTpred)
  par_tab = par_tab[, c("Param", "Param_stan")]
  return(par_tab)
}



