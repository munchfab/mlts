#' Add custom prior distributions to mlts model
#'
#' @param model `data.frame`. Output of \code{\link[mlts]{mlts_model}}.
#' @param default Logical. If set to `TRUE`, default prior specifications are
#' added.
#'
#' @return An object of class `data.frame`.
#' @noRd
#'
mlts_model_priors <- function(model, default = FALSE){

  if(default == TRUE){

    # initialise columns
    cols = c("prior_type", "prior_location", "prior_scale")
    model[,cols] = NA

    # STRUCTURAL MODEL =========================================================
    ## Fixed effects
    model[model$Type=="Fixed effect" & startsWith(model$Param_Label,prefix = "Trait"),cols] = data.frame("normal", 0, 10)
    model[model$Type=="Fixed effect" & model$Param_Label=="Dynamic", cols] = data.frame("normal", 0, 2)
    model[model$Type=="Fixed effect" & model$Param_Label == "Log Innovation Variance",cols] = data.frame("normal", 0, 10)
    model[model$Type=="Fixed effect" & model$Param_Label == "Innovation Variance",cols] = data.frame("cauchy", 0, 2.5)

    # innovation covariance - if random:
    model[model$Type=="Fixed effect" & model$Param_Label == "Log Innovation Covariance",cols] = data.frame("normal",0,10)
    # innovation covariance - if constant:
    model[model$Type=="Fixed effect" & model$Param_Label == "Innovation correlation",cols] = data.frame("LKJ",1,NA)

    ## Random effects
    model[model$Type=="Random effect SD",cols] = data.frame("cauchy", 0, 2.5)
    model[model$Type=="RE correlation",cols] = data.frame("LKJ", 1, NA)

    # MEASUREMENT MODEL
    model[model$Type=="Loading", cols] = data.frame("normal", 1, 0.5)
    model[model$Type=="Item intercepts", cols] = data.frame("normal", 0, 10)
    model[model$Type=="Measurement Error SD", cols] = data.frame("cauchy", 0, 2.5)

    # BETWEEN-LEVEL REGRESSIONS
    model[model$Type=="RE prediction", cols] = data.frame("normal", 0, 10)
    model[model$Type=="Outcome prediction", cols] = data.frame("normal", 0, 10)
    model[model$Type=="Outcome prediction" &
             model$Param_Label == "Residual SD", cols] = data.frame("cauchy", 0, 2.5)
  }


  return(model)

}

