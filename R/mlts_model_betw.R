#' Title
#'
#' @param model data.frame. Output of `model()`-functions.
#' @param ranef_pred character vector or a named list. Include between-level covariate(s)
#' as predictor(s) of all random effects in `model` by entering a vector of unique variable
#' names. Alternatively, to include between-level covariates or differing sets of
#' between-level covariates as predictors of specific random effects, a named
#' list (using the `param`-labels in `model`) can be entered (see details).
#' Note that if a named list is provided, all names that do not match random
#' parameters in `model` will be ignored.
#' @param out_pred character vector or a named list. Include between-level covariate(s)
#' as predictor(s) of all random effects in `model` by entering a vector of unique variable
#' names. Alternatively, to include between-level covariates or differing sets of
#' between-level covariates as predictors of specific random effects, a named
#' list (using the `param`-labels in `model`) can be entered (see details).
#' @param out_pred_add_btw tba
#' @return An object of class `data.frame`.
#' @export
#'
mlts_model_betw <- function(model,ranef_pred = NULL, out_pred=NULL, out_pred_add_btw = NULL){



  ##### ACHTUNG AKTUELL IST ES NOCH MÖGLICH FIXED PARAMETERS IN DIE OUTCOME
  ##### PREDICTION ZU INTEGRIEREN


  # get number of random effects in model
  n_random = nrow(model[model$Type == "Fix effect" & model$isRandom == 1,])
  pars_random = model[model$Type == "Fix effect" & model$isRandom == 1,"Param"]

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
    RE.PRED = mlts_model_priors(RE.PRED, default = T)

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
    OUT.PRED = mlts_model_priors(OUT.PRED, default = T)

    # add to model
    model = dplyr::bind_rows(model, OUT.PRED)
  }

  return(model)
}



