#' Add measurement model structure to a `VARmodel`-object
#'
#' Add (or replace) a measurement model of an existing `VARmodel`-object. As default
#' option, a multiple-indicator model will be included assuming a common between-level factor.
#'
#' @param VARmodel data.frame. Output of VARmodel-Functions.
#' @param q integer. The number of time-varying constructs.
#' @param p integer. A vector of length `q` with the number of manifest
#' indicators per construct.
#' @param btw.factor Logical. If `TRUE` (the default), a common between-level factor
#' is modeled across all indicator variables. If `FALSE`, instead of a between-level
#' factor, indicator mean levels will be included as individual (random) effects stemming
#' from a joint multivariate normal distribution.
#' @param btw.model A list to indicate for which manifest indicator variables a common
#' between-level factor should be modeled (see Details for detailed instructions).
#' At this point restricted to one factor per latent construct.
#' @return An object of class `data.frame`.
#' @details Update a `VARmodel`-object.
#' @examples
#' # build a manifest two-level AR(1) model
#' VARmodel <- VARmodelBuild(q = 1)
#'
#' # add measurement model using three indicators
#' VARmodel <- VARmodelMeasurement(VARmodel, q = 1, p = 3, btw.factor = TRUE)
#'
#' # only indicator-specific (random) intercepts are modeled
#' VARmodel <- VARmodelMeasurement(VARmodel, q = 1, p = 3, btw.factor = TRUE)
#'
#' # A more fine-grained between-level measurement model can be specified via:
#' VARmodel <- VARmodelMeasurement(
#'               VARmodel, q = 1, p = 3,
#'               btw.model = list(""))
#' Which models a common latent factor on the between-level for the first three
#' indicators and a random indicator mean for the fourth indicator.
#' @export
#'
VARmodelMeasurement <- function(VARmodel, q, p, btw.factor = TRUE, btw.model = NULL){

  if(length(p) == 1){
    p = rep(p, times = q)
    if (q > 1) {
      warning("Note: The number of indicators is assumed to be ", p[1],
              " for each latent variable. If this is not intended, please",
              " specify a vector of length q containing the number of",
              " indicators for each latent construct",
              " (see Vignettes for examples).")
    }
  }

  # print a warning if VARmodel already contains a measurement model
  if("Measurement" %in% VARmodel$Model){
    message("VARmodel already contains a measurement model specification which will be overwritten.")
    VARmodel = VARmodel[VARmodel$Model != "Measurement",]
  }


  # print a warning, if model already contains a between-level structure
  if(sum(grepl(VARmodel$Type, pattern = "prediction"))>0){
    warning("Between-level prediction models are removed when changes to the measurement model are made.")
    VARmodel = VARmodel[grepl(VARmodel$Type, pattern = "prediction"),]
  }

  ## build the btw.model-argument if none is provided
  if(is.null(btw.model)){
    btw.model = list()
    for(i in 1:q){
      if(btw.factor == T){
        btw.model[[i]] = c(1:p[i])
      } else {
        btw.model[[i]] = c(NA)
      }
    }
  }

  ## extract information from btw.model in loop over dimensions --------------
  btw.pars = list() # declare a list object to store results
  for(i in 1:q){
    N_inds = p[i]
    N_inds_means = N_inds - length(na.omit((btw.model[[i]])))
    inds_means = which(!(1:p[i] %in% btw.model[[i]]))
    N_etaB_inds = N_inds - N_inds_means
    etaB_inds = which((1:p[i] %in% btw.model[[i]]))
    N_etaB = ifelse(N_inds == N_inds_means, 0, 1)

    # create fixed effect structural part  ==================================
    if(N_etaB == 1){
      etaB.fix = data.frame(
        "Model" = "Structural",
        "Level" = "Within",
        "Type" = "Fix effect",
        "Param" = paste0("etaB_",i),
        "Param_Label" = "Trait (latent factor)",
        "isRandom" = 1,
        "Constraint"= c(NA)
      )
      etaB.rand = data.frame(
        "Model" = "Structural",
        "Level" = "Between",
        "Type" = "Random effect SD",
        "Param" = paste0("sigma_etaB_",i),
        "Param_Label" = "Trait (latent factor)",
        "isRandom" = 0,
        "Constraint"= c(NA)
      )
    }
    if(N_inds_means > 0){
      mus.fix = data.frame(
        "Model" = "Structural",
        "Level" = "Within",
        "Type" = "Fix effect",
        "Param" = paste0("mu_",i,".",inds_means),
        "Param_Label" = "Trait (indicator mean)",
        "isRandom" = 1,
        "Constraint"= c(NA)
      )
      mus.rand = data.frame(
        "Model" = "Structural",
        "Level" = "Between",
        "Type" = "Random effect SD",
        "Param" = paste0("sigma_mu_",i,".",inds_means),
        "Param_Label" = "Trait (indicator mean)",
        "isRandom" = 0,
        "Constraint"= c(NA)
      )
    }

    # combine and replace in VARmodel
    if(N_etaB == 1 & N_inds_means > 0){
      fix = rbind(etaB.fix, mus.fix)
      rand = rbind(etaB.rand, mus.rand)
    } else if(N_etaB == 1 & N_inds_means == 0){
      fix = etaB.fix
      rand = etaB.rand
    } else if(N_etaB == 0 & N_inds_means > 0){
      fix = mus.fix
      rand = mus.rand
    }
    # replacement
    fix = VARmodelPriors(fix, default = T)    # add default priors
    rand = VARmodelPriors(rand, default = T)

    row_to_repl = which(VARmodel$Param == paste0("mu_",i) & VARmodel$Type == "Fix effect")
    VARmodel = replace_model_row(VARmodel, row_to_repl, fix)
    row_to_repl = which(VARmodel$Param == paste0("sigma_mu_",i) & VARmodel$Type == "Random effect SD")
    VARmodel = replace_model_row(VARmodel, row_to_repl, rand)

    # update the random effect correlations
    VARmodel = update_model_REcors(VARmodel)


    # Between-level =============================================================
    ## create additional rows
    if(N_etaB != 0){
      ## Item-Intercepts (first fixed to 1 as default option)
      itemInts = data.frame(
        "Model" = "Measurement",
        "Level" = "Between",
        "Type" = "Item intercepts",
        "Param"= c(paste0("alpha_",i,".",etaB_inds)),
        "Param_Label" = c(""),
        "Constraint" = c("= 0", rep("free", (N_etaB_inds-1)))
      )

      ## Loadings (first fixed to 1 as default option)
      loadsB = data.frame(
        "Model" = "Measurement",
        "Level" = "Between",
        "Type" = "Loading",
        "Param"= c(paste0("lambdaB_",i,".",etaB_inds)),
        "Param_Label" = c(""),
        "Constraint" = c("= 1", rep("free", (N_etaB_inds-1))))

      ## Residual variances
      errVarB = data.frame(
        "Model" = "Measurement",
        "Level" = "Between",
        "Type" = "Measurement Error SD",
        "Param"= c(paste0("sigmaB_",i,".",etaB_inds)),
        "Param_Label" = c(""),
        "Constraint" = c("= 0", rep("free", (N_etaB_inds-1)))
      )

      # combine
      btw.pars[[i]] = rbind(itemInts, loadsB, errVarB)
    } else {
      btw.pars[[i]] = NULL
    }
  }
  # combine
  between = as.data.frame(data.table::rbindlist(btw.pars))



  # Within-level =============================================================
  ## create additional rows
  # for each of the q constructs depending on the number of indicators:
  mm.pars = list()
  for(i in 1:q){
    ## Loadings (first fixed to 1 as default option)
    loadsW = data.frame(
      "Model" = "Measurement",
      "Level" = "Within",
      "Type" = "Loading",
      "Param"= c(paste0("lambdaW_",i,".",1:p[i])),
      "Param_Label" = c(""),
      "Constraint" = c("= 1", rep("free", (p[i]-1)))
    )

    ## Residual variances
    errVarW = data.frame(
      "Model" = "Measurement",
      "Level" = "Within",
      "Type" = "Measurement Error SD",
      "Param"= c(paste0("sigmaW_",i,".",1:p[i])),
      "Param_Label" = c(""),
      "Constraint" = ifelse(p[i] == 1, "= 0", "free")
    )
    mm.pars[[i]] = plyr::rbind.fill(loadsW, errVarW)
  }

  # combine
  within = as.data.frame(data.table::rbindlist(mm.pars))
  mm.pars = rbind(between, within)
  mm.pars = dplyr::arrange(mm.pars, Level, Type, Param)

  ## add default priors
  mm.pars = VARmodelPriors(VARmodel = mm.pars, default = T)


  # add to structural part
  VARmodel = plyr::rbind.fill(VARmodel, mm.pars)



  return(VARmodel)

}



