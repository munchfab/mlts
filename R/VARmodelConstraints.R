#' Title
#'
#' @param VARmodel data.frame. Output of VARmodel-Functions.
#' @param FiVARmodeledInnoVars logical. FiVARmodel all random effect variances (eVARmodelcept those
#' of individual traits) to zero.
#' @param FiVARmodeledCovs logical. Set all innovation covariances to a constant value.
#' @param CovsZero logical. Set to TRUE to treat all innovations as independent.
#' @param FEis0 character. A character vector to indeVARmodel which fiVARmodeled model parameters
#' should be fiVARmodeled to zero (Note: this results in removing the random effect
#' variance of the respective parameter).
#' @param REis0 logical. Set to TRUE to treat all innovations as independent.
#'
#' @return An object of class `data.frame`.
#' @export
#'
VARmodelConstraints <- function(VARmodel, FixDynamics = F, FixInnoVars = F,
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




  # update random effect correlations
  rand.pars = (VARmodel[VARmodel$Type == "Fix effect" & VARmodel$isRandom == 1,"Param"])
  n_rand = length(rand.pars)
  btw.cov_pars = c()
  if(n_rand>1){
    n_cors = (n_rand * (n_rand-1))/2
    qs = c()
    ps = c()
    for(i in 1:(n_rand-1)){
    qs = c(qs, rep(rand.pars[i], each = n_rand-i))
    }
    for(i in 2:n_rand){
      ps = c(ps, rep(rand.pars[i:n_rand], 1))
    }

  btw.cov_pars = paste0("r_", qs,".", ps)

  # remove and replace
  VARmodel = VARmodel[VARmodel$Type!="RE correlation" |
    (VARmodel$Type=="RE correlation" & VARmodel$Param %in% btw.cov_pars),]
  } else if(n_rand == 1){

    VARmodel = VARmodel[VARmodel$Type != "RE correlation",]

  }


  # update priors =============================================================
  VARmodel = VARmodelPriors(VARmodel = VARmodel, default = T)


  return(VARmodel)

}
