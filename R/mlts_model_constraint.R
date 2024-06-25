#' Add parameter constraints to mlts model
#'
#' @param model `data.frame`. Output of \code{\link[mlts]{mlts_model}}.
#' @param fix_dynamics Logical. Fix all random effect variances of autoregressive and
#' cross-lagged effects to zero (constraining parameters to be equal across clusters).
#' @param fix_inno_vars Logical. Fix all random effect variances of innovation variances
#' to zero (constraining parameters to be equal across clusters).
#' @param fix_inno_covs Logical. Fix all random effect variances of innovation covariances
#' to zero (constraining parameters to be equal across clusters).
#' @param inno_covs_zero Logical. Set to `TRUE` to treat all innovations as independent.
#' @param fixef_zero Character. A character vector to index which fixed effects
#' (referring to the parameter labels in `model$Param`) should be constrained to zero
#' (Note: this also results in removing the random effect variance of the respective parameter).
#' @param ranef_zero Character. A character vector to index which random effect variances
#' (referring to the parameter labels in `model$Param`) should be constrained to zero.
#'
#' @return An object of class `data.frame`.
#' @noRd
#'
#' @examples
#' # simple vector-autoregressive mlts model with q = 2 time-series constructs
#' # and autoregressive max lag of second order
#' var_model <- mlts_model(q = 2, max_lag = 2)
#'
#' # to fix only specific parameters, provide a vector with names
#' # referring to model parameters
#' # e.g., fix the first-order autoregressive effect of the second
#' # construct (by setting its variance to zero)
#' var_model1 <- mlts_model_constraint(
#'   model = var_model, ranef_zero = "phi(1)_22"
#' )
#' # Parameters fixed to zero (i.e., the variance of phi(1)_22) are
#' # deleted from the model data frame
#'
#' # innovation variances can also be set to zero
#' # e.g., fix first-order autoregressive effect and log innovation
#' # variance of first construct
#' var_model2 <- mlts_model_constraint(
#'   model = var_model,
#'   ranef_zero = c("phi(1)_22", "ln.sigma2_1")
#' )
#'
#' # if all autoregressive and cross-lagged parameters should be fixed,
#' # fix_dynamics = TRUE can be used
#' var_model3 <- mlts_model_constraint(model = var_model, fix_dynamics = TRUE)
#'
#' # if all innovation variances should be fixed,
#' # fix_inno_vars = TRUE can be used
#' var_model3 <- mlts_model_constraint(model = var_model, fix_inno_vars = TRUE)
#'
#'
#' # fixed effects can be set to zero by providing a vector with
#' # parameter names to fixef_zero
#' # e.g., fixing the cross-lagged effect of first order from the first
#' # construct to the second construct to zero
#' var_model5 <- mlts_model_constraint(
#'   model = var_model, fixef_zero = "phi(1)_21"
#' )
#'

mlts_model_constraint <- function(model, fix_dynamics = FALSE, fix_inno_vars = FALSE,
                                  fix_inno_covs = FALSE, inno_covs_zero = FALSE,
                                  fixef_zero = NULL, ranef_zero = NULL
){

  if(fix_dynamics == TRUE){

    # update indentifier column
    model$isRandom[grepl(model$Param_Label, pattern = "Dynamic")] = 0

    # remove random effect SDs
    model = model[!(model$Type == "Random effect SD" &
                            model$Param_Label == "Dynamic"),]
  }


  if(fix_inno_vars == TRUE){

    # update indentifier column
    model$isRandom[grepl(model$Param, pattern = "ln.sigma2_")] = 0

    # remove random effect SDs
    model = model[!(model$Type == "Random effect SD" &
                            model$Param_Label == "Log Innovation Variance"),]
    # adjust labels of fixed effects
    model[(model$Type == "Fixed effect" &
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
    model$Param_Label[(model$Type == "Fixed effect" &
                            model$Param_Label == "Log Innovation Covariance")] = "Innovation correlation"

    # adjust params
    model$Param[grepl(model$Param, pattern = "ln.sigma_")] = gsub(
      model$Param[grepl(model$Param, pattern = "ln.sigma_")],
      pattern = "ln.sigma_", replacement = "r.zeta_")
  }

  if(inno_covs_zero == TRUE){
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
    model$isRandom[model$Type=="Fixed effect" & model$Param %in% c(ranef_zero)] = 0

    isFixed = paste0("sigma_", model$Param[model$Type == "Fixed effect" & model$isRandom == 0])
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
  model = mlts_model_priors(model = model, default = TRUE)


  return(model)

}
