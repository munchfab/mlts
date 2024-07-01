#' Add between-level variables to mlts model
#'
#' @param model `data.frame`. Output of \code{\link[mlts]{mlts_model}}.
#' @param ranef_pred A character vector or a named list. Include between-level covariate(s)
#' as predictor(s) of all random effects in `model` by entering a vector of unique variable
#' names. Alternatively, to include between-level covariates or differing sets of
#' between-level covariates as predictors of specific random effects, a named
#' list (using the labels in `model$Param`) can be entered (see examples).
#' Note that if a named list is provided, all names that do not match random
#' parameters in `model` will be ignored.
#' @param out_pred A character vector or a named list. Include between-level covariate(s)
#' as predictor(s) of all random effects in `model` by entering a vector of unique variable
#' names. Alternatively, to include between-level covariates or differing sets of
#' between-level covariates as predictors of specific random effects, a named
#' list (using the labels in `model$Param`) can be entered (see examples). Note that
#' if a named list is provided, all character strings in the vector of each list
#' (with independent variables) element that do not match random effect parameter names
#' in `model$Param` will be treated as additional between-level predictors.
#' @param out_pred_add_btw A character vector. All inputs will be treated as
#' between-level covariates to be used as additional predictors of all outcomes specified
#' in `out_pred`.
#' @return An object of class `data.frame`.
#' @noRd
#'
#' @examples
#' # simple autoregressive mlts model with q = 2 time-series constructs
#' ar_model <- mlts_model(q = 2)
#'
#' # add between-level variables to predict all random parameters in the model
#' ar_model1 <- mlts_model_betw(model = ar_model, ranef_pred = "x")
#'
#' # if only some of the random parameters should be predicted,
#' # a named list can be provided of the format `list("parameter" = "covariate")`
#' # e.g., to only predict the random intercept (trait) of the first construct
#' ar_model2 <- mlts_model_betw(
#'   model = ar_model, ranef_pred = list("mu_1" = "x")
#' )
#'
#' # predict outcomes using (all) random parameters in the model
#' ar_model3 <- mlts_model_betw(model = ar_model, out_pred = "y")
#'
#' # to use only specific random parameter, provide a named list
#' # e.g., predict y by random intercept (trait) and autoregressive effect
#' # of first construct
#' ar_model4 <- mlts_model_betw(
#'   model = ar_model, out_pred = list("y" = c("mu_1", "phi(1)_11"))
#' )
#'
#'
mlts_model_betw <- function(model,ranef_pred = NULL, out_pred=NULL, out_pred_add_btw = NULL){

  # print a warning, if specific outcome prediction model is introduced,
  # and additional between-level covariates shoueld be entered:
  if(is.list(out_pred) & !is.null(out_pred_add_btw)){
    warning(
      "If a list input is provided for `out_pred`, additional between-level
      covariates to be used as outcome predictors entered via `out_pred_add_btw`
      will be ignored. You can add those variables to the list of named
      vectors in `out_pred`."
      )
  }


  # get number of random effects in model
  n_random = nrow(model[model$Type == "Fixed effect" & model$isRandom == 1,])
  pars_random = model[model$Type == "Fixed effect" & model$isRandom == 1,"Param"]

  # offer different input options:
  ## Wenn input einzelner character oder vector mit charactern ist, verwende
  ## Kovariate als Prädiktor für alle REs, ansonsten named list um spezifische
  ## Modelle mit variierenden Prädiktoren für einzelne REs zu definieren

  # time-invariant covariates as predictors of random effects
  if(is.list(ranef_pred) == FALSE){
    # use covariates as predictors of all random effects as default
    b_pred_params = c()
    for(j in 1:n_random){
      for(i in 1:length(ranef_pred)){
        # create parameter names
        b_pred_params = c(b_pred_params, paste0("b_",pars_random[j], ".ON.",ranef_pred[i]))
      }
    }
  } else {
    b_pred_params = c()
    n_rand_preds = length(ranef_pred)
    rand_preds_params = pars_random[pars_random %in% names(ranef_pred)]
    for(j in 1:n_rand_preds){
      # get matching position on ranef_pred list
      k = which(names(ranef_pred) == rand_preds_params[j])
      for(i in 1:length(ranef_pred[[k]])){
        # create parameter names
        b_pred_params = c(b_pred_params, paste0("b_",rand_preds_params[j], ".ON.",ranef_pred[[k]][i]))
      }
    }
  }

  # combine a information in a data frame
  if(!is.null(ranef_pred)){
    RE.PRED = data.frame(
      "Model" = "Structural",
      "Level" = "Between",
      "Type" = "RE prediction",
      "Param" = b_pred_params,
      "Param_Label" = "regression weight",
      "isRandom" = 0
    )

    # add priors
    RE.PRED = mlts_model_priors(RE.PRED, default = TRUE)

    # add to model
    model = dplyr::bind_rows(model, RE.PRED)
  }

  # outcome prediction models
  if(is.list(out_pred) == FALSE){
    # use all random effects as predictors as default option
    if(!is.null(out_pred_add_btw)){
      # add additional between-level variables to predictors
      n_random = n_random + length(out_pred_add_btw)
      pars_random = c(pars_random, out_pred_add_btw)
    }

    out_pred_params = c()
    for(i in 1:length(out_pred)){
      for(j in 1:n_random){
        # create parameter names of regression coefficients
        out_pred_params = c(out_pred_params, paste0("b_",out_pred[i], ".ON.",pars_random[j]))
      }
      # add intercept and SD of outcome residual variance
      out_pred_params = c(
        paste0("alpha_",out_pred[i]),
        out_pred_params,
        paste0("sigma_",out_pred[i]))
    }
  } else {
    out_pred_params = c()
    n_out_pred = length(out_pred)
    for(j in 1:n_out_pred){
      # reorder entered REs as predictors in order of random effects
      out_pred[[j]] =  c(pars_random[pars_random %in% out_pred[[j]]], #REs ordered
                         out_pred[[j]][!(out_pred[[j]] %in% pars_random)])

      for(i in 1:length(out_pred[[j]])){
        # create parameter names of regression coefficients
        out_pred_params = c(out_pred_params, paste0("b_",names(out_pred)[j], ".ON.",out_pred[[j]][i]))
      }
      # add intercept and SD of outcome residual variance
      out_pred_params = c(paste0("alpha_",names(out_pred)[j]),
                          out_pred_params,
                          paste0("sigma_",names(out_pred)[j]))
    }
  }

  # combine a information in a data frame
  if(!is.null(out_pred)){
    OUT.PRED = data.frame(
      "Model" = "Structural",
      "Level" = "Between",
      "Type" = "Outcome prediction",
      "Param" = out_pred_params,
      "Param_Label" = ifelse(startsWith(out_pred_params,"sigma_"),
                             "Residual SD", "regression weight"),
      "isRandom" = 0
    )
    # intercept
    OUT.PRED$Param_Label = ifelse(
      startsWith(OUT.PRED$Param, "alpha_"),"intercept",OUT.PRED$Param_Label)

    # add priors
    OUT.PRED = mlts_model_priors(OUT.PRED, default = TRUE)

    # add to model
    model = dplyr::bind_rows(model, OUT.PRED)
  }

  return(model)
}



