#' Add parameter constraints to mlts model
#'
#' @param model data.frame. Output of \code{\link[mlts]{mlts_model}}.
#' @param fix_dynamics logical. Fix all random effect variances (except those
#' of individual traits) to zero.
#' @param fix_inno_vars logical. Set all innovation variances to a constant value.
#' @param fix_inno_covs logical. Set all innovation covariances to a constant value.
#' @param inno_covs_zero logical. Set to TRUE to treat all innovations as independent.
#' @param fixef_zero character. A character vector to index which fixed model parameters
#' should be fixed to zero (Note: this results in removing the random effect
#' variance of the respective parameter).
#' @param ranef_zero logical. Set to TRUE to treat all innovations as independent.
#'
#' @return An object of class `data.frame`.
#' @export
#'
mlts_model_constraint <- function(model, fix_dynamics = F, fix_inno_vars = F,
                                fix_inno_covs = F, inno_covs_zero = F,
                                fixef_zero = NULL, ranef_zero = NULL
){

  if(fix_dynamics == T){

    # update indentifier column
    model$isRandom[grepl(model$Param_Label, pattern = "Dynamic")] = 0

    # remove random effect SDs
    model = model[!(model$Type == "Random effect SD" &
                            model$Param_Label == "Dynamic"),]
  }


  if(fix_inno_vars == T){

    # update indentifier column
    model$isRandom[grepl(model$Param, pattern = "ln.sigma2_")] = 0

    # remove random effect SDs
    model = model[!(model$Type == "Random effect SD" &
                            model$Param_Label == "Log Innovation Variance"),]
    # adjust labels of fixed effects
    model[(model$Type == "Fix effect" &
                model$Param_Label == "Log Innovation Variance"), "Param_Label"] = "Innovation Variance"
    # adjust params
    model$Param[grepl(model$Param, pattern = "ln.sigma2_")] = gsub(
      model$Param[grepl(model$Param, pattern = "ln.sigma2_")],
      pattern = "ln.sigma2_", replacement = "sigma_")

  }

  if(fix_inno_covs == T){

    # update indentifier column
    model$isRandom[model$Param_Label == "Log Innovation Covariance"] = 0

    # remove random effect SDs
    model = model[!(model$Type == "Random effect SD" &
                            model$Param_Label == "Log Innovation Covariance"),]

    # adjust labels of fixed effects
    model$Param_Label[(model$Type == "Fix effect" &
                            model$Param_Label == "Log Innovation Covariance")] = "Innovation correlation"

    # adjust params
    model$Param[grepl(model$Param, pattern = "ln.sigma_")] = gsub(
      model$Param[grepl(model$Param, pattern = "ln.sigma_")],
      pattern = "ln.sigma_", replacement = "r.zeta_")
  }

  if(inno_covs_zero == T){
    model = model[!grepl(model$Param_Label, pattern = "Covariance"),]
    model = model[!grepl(model$Param, pattern = "r.zeta"),]
  }

  # remove individual effects from the dynamic model
  if(!is.null(fixef_zero)){
    model = model[!(model$Param %in% c(fixef_zero)),]
    # remove respective random effects
    model = model[!(model$Param %in% c(paste0("sigma_",fixef_zero))),]
  }

  # remove random effect SDs of constant parameters
  if(!is.null(ranef_zero)){
    # update identifier column
    model$isRandom[model$Type=="Fix effect" & model$Param %in% c(ranef_zero)] = 0

    isFixed = paste0("sigma_", model$Param[model$Type == "Fix effect" & model$isRandom == 0])
    model = model[!(model$Type == "Random effect SD" & model$Param %in% isFixed),]

    # adjust parameter labels for constant innovation variances
    model[model$Param %in% ranef_zero, "Param_Label"] = gsub(
      model[model$Param %in% ranef_zero, "Param_Label"],
      pattern = "Log Innovation Variance", replacement = "Innovation Variance")

    # adjust parameter labels for constant innovation covariances
    model[model$Param %in% ranef_zero, "Param_Label"] = gsub(
      model[model$Param %in% ranef_zero, "Param_Label"],
      pattern = "Log Innovation Covariance", replacement = "Innovation correlation")

    # adjust params
    model[model$Param %in% ranef_zero, "Param"] = gsub(
      model[model$Param %in% ranef_zero, "Param"],
      pattern = "ln.sigma2_", replacement = "sigma_")
    model[model$Param %in% ranef_zero, "Param"] = gsub(
      model[model$Param %in% ranef_zero, "Param"],
      pattern = "ln.sigma_", replacement = "r.zeta_")
  }

  # update RE correlations
  model = update_model_REcors(model)

  # update priors =============================================================
  model = mlts_model_priors(model = model, default = T)


  return(model)

}
