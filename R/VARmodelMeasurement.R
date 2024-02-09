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
#' across is modeled across all indicator variables. If `FALSE`, instead of a between-level
#' factor, indicator mean levels will be included as individual (random) effects stemming
#' from a joint multivariate normal distribution.
#' @param btw.model A list to indicate for which manifest indicator variables a common
#' between-level factor should be modeled (see Details for detailed instructions). At this point restricted to one factor per latent construct.
#' @return An object of class `data.frame`.
#' @details Update a `VARmodel`-object.
#' @examples
#' # build a manifest two-level AR(1) model
#' VARmodel <- VARmodelBuild(q = 1)
#'
#' # add measurement model using three indicators
#' VARmodel <- VARmodelMeasurement(VARmodel, q = 1, p = 3, btw.factor = TRUE)
#' @export
#'
VARmodelMeasurement <- function(VARmodel, q, p, btw.factor = TRUE, btw.model = NULL){


  if(length(p) == 1){
    p = rep(p, times = q)
  }

  # print a warning if VARmodel already contains a measurement model
  if("Measurement" %in% VARmodel$Model){
    message("VARmodel already contains a measurement model specification which will be overwritten.")
    VARmodel = VARmodel[VARmodel$Model != "Measurement",]
    }





  # for each of the q constructs depending on the number of indicators:
  mm.pars = list()
  for(i in 1:q){

    # Within-level
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

    # Between-level
    ## Item-Intercepts (first fixed to 1 as default option)
    itemInts = data.frame(
      "Model" = "Measurement",
      "Level" = "Between",
      "Type" = "Item intercepts",
      "Param"= c(paste0("alpha_",i,".",1:p[i])),
      "Param_Label" = c(""),
      "Constraint" = c("= 0", rep("free", (p[i]-1)))
    )

    ## Loadings (first fixed to 1 as default option)
    loadsB = data.frame(
      "Model" = "Measurement",
      "Level" = "Between",
      "Type" = "Loading",
      "Param"= c(paste0("lambdaB_",i,".",1:p[i])),
      "Param_Label" = c(""),
      "Constraint" = c("= 1", rep("free", (p[i]-1)))
    )

    ## Residual variances
    errVarB = data.frame(
      "Model" = "Measurement",
      "Level" = "Between",
      "Type" = "Measurement Error SD",
      "Param"= c(paste0("sigmaB_",i,".",1:p[i])),
      "Param_Label" = c(""),
      "Constraint" = c("= 0", rep("free", (p[i]-1)))
    )

    mm.pars[[i]] = rbind(
      loadsW, errVarW, itemInts, loadsB, errVarB
    )

    }

  # combine
  mm.pars = as.data.frame(data.table::rbindlist(mm.pars))
  mm.pars = dplyr::arrange(mm.pars,Level, Type, Param)

  ## add default priors
  mm.pars = VARmodelPriors(VARmodel = mm.pars, default = T)


  # add to structural part
  VARmodel = plyr::rbind.fill(VARmodel, mm.pars)


  return(VARmodel)

}
