#' Title
#'
#' @param VARmodel data.frame. Output of VARmodel-Functions.
#' @param default logical. If set to `TRUE`, default prior specifications are 
#' added.
#'
#' @return An object of class `data.frame`.
#' @export
VARmodelPriors <- function(VARmodel, default = F){
  
  if(default == T){
    
    # initialise columns 
    cols = c("prior_type", "prior_location", "prior_scale")
    VARmodel[,cols] = NA
    
    # STRUCTURAL MODEL =========================================================
    ## Fixed effects 
    VARmodel[VARmodel$Type=="Fix effect" & VARmodel$Param_Label=="Trait",cols] = data.frame("normal", 0, 10)
    VARmodel[VARmodel$Type=="Fix effect" & VARmodel$Param_Label=="Dynamic", cols] = data.frame("normal", 0, 2)
    VARmodel[VARmodel$Type=="Fix effect" & VARmodel$Param_Label == "Log Innovation Variance",cols] = data.frame("normal", 0, 10)
    VARmodel[VARmodel$Type=="Fix effect" & VARmodel$Param_Label == "Innovation Variance",cols] = data.frame("cauchy", 0, 2.5)
    
    # innovation covariance - if random: 
    VARmodel[VARmodel$Type=="Fix effect" & VARmodel$Param_Label == "Log Innovation Covariance",cols] = data.frame("normal",0,10)
    # innovation covariance - if constant: 
    VARmodel[VARmodel$Type=="Fix effect" & VARmodel$Param_Label == "Innovation correlation",cols] = data.frame("LKJ",1,NA)

    ## Random effects 
    VARmodel[VARmodel$Type=="Random effect SD",cols] = data.frame("cauchy", 0, 2.5)
    VARmodel[VARmodel$Type=="Random effect correlation",cols] = data.frame("LKJ", 1, NA)
    
    # MEASUREMENT MODEL 
    VARmodel[VARmodel$Type=="Loading", cols] = data.frame("normal", 1, 0.5)
    VARmodel[VARmodel$Type=="Item intercepts", cols] = data.frame("normal", 0, 10)
    VARmodel[VARmodel$Type=="Measurement Error SD", cols] = data.frame("cauchy", 0, 2.5)
  }
  
  
  return(VARmodel)
  
}

