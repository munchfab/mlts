#' Title
#'
#' @param VARmodel data.frame. Output of VARmodel-Functions.
#' @param FixDynamics logical. Fix all random effect variances (except those
#' of individual traits) to zero.
#' @param FixInnoVars logical. Set all innovation variances to a constant value.
#' @param FixInnoCovs logical. Set all innovation covariances to a constant value.
#' @param InnoCovsZero logical. Set to TRUE to treat all innovations as independent.
#' @param FEis0 character. A character vector to indeVARmodel which fiVARmodeled model parameters
#' should be fiVARmodeled to zero (Note: this results in removing the random effect
#' variance of the respective parameter).
#' @param REis0 logical. Set to TRUE to treat all innovations as independent.
#'
#' @return An object of class `data.frame`.
#' @export
#'
mlts_model_constraint <- function(VARmodel, FixDynamics = F, FixInnoVars = F,
                                FixInnoCovs = F, InnoCovsZero = F,
                                FEis0 = NULL, REis0 = NULL
){

  if(FixDynamics == T){

    # update indentifier column
    VARmodel$isRandom[grepl(VARmodel$Param_Label, pattern = "Dynamic")] = 0

    # remove random effect SDs
    VARmodel = VARmodel[!(VARmodel$Type == "Random effect SD" &
                            VARmodel$Param_Label == "Dynamic"),]
  }


  if(FixInnoVars == T){

    # update indentifier column
    VARmodel$isRandom[grepl(VARmodel$Param, pattern = "ln.sigma2_")] = 0

    # remove random effect SDs
    VARmodel = VARmodel[!(VARmodel$Type == "Random effect SD" &
                            VARmodel$Param_Label == "Log Innovation Variance"),]
    # adjust labels of fixed effects
    VARmodel[(VARmodel$Type == "Fix effect" &
                VARmodel$Param_Label == "Log Innovation Variance"), "Param_Label"] = "Innovation Variance"
    # adjust params
    VARmodel$Param[grepl(VARmodel$Param, pattern = "ln.sigma2_")] = gsub(
      VARmodel$Param[grepl(VARmodel$Param, pattern = "ln.sigma2_")],
      pattern = "ln.sigma2_", replacement = "sigma_")

  }

  if(FixInnoCovs == T){

    # update indentifier column
    VARmodel$isRandom[VARmodel$Param_Label == "Log Innovation Covariance"] = 0

    # remove random effect SDs
    VARmodel = VARmodel[!(VARmodel$Type == "Random effect SD" &
                            VARmodel$Param_Label == "Log Innovation Covariance"),]

    # adjust labels of fixed effects
    VARmodel$Param_Label[(VARmodel$Type == "Fix effect" &
                            VARmodel$Param_Label == "Log Innovation Covariance")] = "Innovation correlation"

    # adjust params
    VARmodel$Param[grepl(VARmodel$Param, pattern = "ln.sigma_")] = gsub(
      VARmodel$Param[grepl(VARmodel$Param, pattern = "ln.sigma_")],
      pattern = "ln.sigma_", replacement = "r.zeta_")
  }

  if(InnoCovsZero == T){
    VARmodel = VARmodel[!grepl(VARmodel$Param_Label, pattern = "Covariance"),]
    VARmodel = VARmodel[!grepl(VARmodel$Param, pattern = "r.zeta"),]
  }

  # remove individual effects from the dynamic model
  if(!is.null(FEis0)){
    VARmodel = VARmodel[!(VARmodel$Param %in% c(FEis0)),]
    # remove respective random effects
    VARmodel = VARmodel[!(VARmodel$Param %in% c(paste0("sigma_",FEis0))),]
  }

  # remove random effect SDs of constant parameters
  if(!is.null(REis0)){
    # update identifier column
    VARmodel$isRandom[VARmodel$Type=="Fix effect" & VARmodel$Param %in% c(REis0)] = 0

    isFixed = paste0("sigma_", VARmodel$Param[VARmodel$Type == "Fix effect" & VARmodel$isRandom == 0])
    VARmodel = VARmodel[!(VARmodel$Type == "Random effect SD" & VARmodel$Param %in% isFixed),]

    # adjust parameter labels for constant innovation variances
    VARmodel[VARmodel$Param %in% REis0, "Param_Label"] = gsub(
      VARmodel[VARmodel$Param %in% REis0, "Param_Label"],
      pattern = "Log Innovation Variance", replacement = "Innovation Variance")

    # adjust parameter labels for constant innovation covariances
    VARmodel[VARmodel$Param %in% REis0, "Param_Label"] = gsub(
      VARmodel[VARmodel$Param %in% REis0, "Param_Label"],
      pattern = "Log Innovation Covariance", replacement = "Innovation correlation")

    # adjust params
    VARmodel[VARmodel$Param %in% REis0, "Param"] = gsub(
      VARmodel[VARmodel$Param %in% REis0, "Param"],
      pattern = "ln.sigma2_", replacement = "sigma_")
    VARmodel[VARmodel$Param %in% REis0, "Param"] = gsub(
      VARmodel[VARmodel$Param %in% REis0, "Param"],
      pattern = "ln.sigma_", replacement = "r.zeta_")
  }

  # update RE correlations
  VARmodel = update_model_REcors(VARmodel)

  # update priors =============================================================
  VARmodel = mlts_model_priors(VARmodel = VARmodel, default = T)


  return(VARmodel)

}
