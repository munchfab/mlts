#' Title
#'
#' @param VARmodel data.frame. Output of `VARmodel()`-functions.
#' @param RE.pred character vector or a named list. Include between-level covariate(s)
#' as predictor(s) of all random effects in `VARmodel` by entering a vector of unique variable
#' names. Alternatively, to include between-level covariates or differing sets of
#' between-level covariates as predictors of specific random effects, a named
#' list (using the `param`-labels in `VARmodel`) can be entered (see details).
#' Note that if a named list is provided, all names that do not match random
#' parameters in `VARmodel` will be ignored.
#' @param out.pred character vector or a named list. Include between-level covariate(s)
#' as predictor(s) of all random effects in `VARmodel` by entering a vector of unique variable
#' names. Alternatively, to include between-level covariates or differing sets of
#' between-level covariates as predictors of specific random effects, a named
#' list (using the `param`-labels in `VARmodel`) can be entered (see details).
#' @return An object of class `data.frame`.
#' @export
#'
VARmodelBetween <- function(VARmodel,RE.pred = NULL, out.pred=NULL, out.pred.add.btw = NULL){



  ##### ACHTUNG AKTUELL IST ES NOCH MÖGLICH FIXED PARAMETERS IN DIE OUTCOME
  ##### PREDICTION ZU INTEGRIEREN


  # get number of random effects in VARmodel
  n_random = nrow(VARmodel[VARmodel$Type == "Fix effect" & VARmodel$isRandom == 1,])
  pars_random = VARmodel[VARmodel$Type == "Fix effect" & VARmodel$isRandom == 1,"Param"]

  # offer different input options:
  ## Wenn input einzelner character oder vector mit charactern ist, verwende
  ## Kovariate als Prädiktor für alle REs, ansonsten named list um spezifische
  ## Modelle mit variierenden Prädiktoren für einzelne REs zu definieren

  # time-invariant covariates as predictors of random effects
  if(is.list(RE.pred) == FALSE){
    # use covariates as predictors of all random effects as default
    b_pred_params = c()
    for(j in 1:n_random){
      for(i in 1:length(RE.pred)){
        # create parameter names
        b_pred_params = c(b_pred_params, paste0("b_",pars_random[j], ".ON.",RE.pred[i]))
      }
    }
  } else {
    b_pred_params = c()
    n_rand_preds = length(RE.pred)
    rand_preds_params = pars_random[pars_random %in% names(RE.pred)]
    for(j in 1:n_rand_preds){
      # get matching position on RE.pred list
      k = which(names(RE.pred) == rand_preds_params[j])
      for(i in 1:length(RE.pred[[k]])){
        # create parameter names
        b_pred_params = c(b_pred_params, paste0("b_",rand_preds_params[j], ".ON.",RE.pred[[k]][i]))
      }
    }
  }

  # combine a information in a data frame
  if(!is.null(RE.pred)){
    RE.PRED = data.frame(
      "Model" = "Structural",
      "Level" = "Between",
      "Type" = "RE prediction",
      "Param" = b_pred_params,
      "Param_Label" = "regression weight",
      "isRandom" = 0
    )

    # add priors
    RE.PRED = VARmodelPriors(RE.PRED, default = T)

    # add to VARmodel
    VARmodel = rbind(VARmodel, RE.PRED)
  }

  # outcome prediction models
  if(is.list(out.pred) == FALSE){
    # use all random effects as predictors as default option
    if(!is.null(out.pred.add.btw)){
      # add additional between-level variables to predictors
      n_random = n_random + length(out.pred.add.btw)
      pars_random = c(pars_random, out.pred.add.btw)
    }

    out_pred_params = c()
    for(i in 1:length(out.pred)){
      for(j in 1:n_random){
        # create parameter names of regression coefficients
        out_pred_params = c(out_pred_params, paste0("b_",out.pred[i], ".ON.",pars_random[j]))
      }
      # add intercept and SD of outcome residual variance
      out_pred_params = c(
        paste0("alpha_",out.pred[i]),
        out_pred_params,
        paste0("sigma_",out.pred[i]))
    }
  } else {
    out_pred_params = c()
    n_out_pred = length(out.pred)
    for(j in 1:n_out_pred){
      # reorder entered REs as predictors in order of random effects
      out.pred[[j]] =  c(pars_random[pars_random %in% out.pred[[j]]], #REs ordered
                         out.pred[[j]][!(out.pred[[j]] %in% pars_random)])

      for(i in 1:length(out.pred[[j]])){
        # create parameter names of regression coefficients
        out_pred_params = c(out_pred_params, paste0("b_",names(out.pred)[j], ".ON.",out.pred[[j]][i]))
      }
      # add intercept and SD of outcome residual variance
      out_pred_params = c(paste0("alpha_",names(out.pred)[j]),
                          out_pred_params,
                          paste0("sigma_",names(out.pred)[j]))
    }
  }

  # combine a information in a data frame
  if(!is.null(out.pred)){
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
    OUT.PRED = VARmodelPriors(OUT.PRED, default = T)

    # add to VARmodel
    VARmodel = rbind(VARmodel, OUT.PRED)
  }

  ###### Should we update the label of random effects SDs if
  ###### RE predictors are provided?

  return(VARmodel)
}



