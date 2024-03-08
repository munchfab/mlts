#' Title
#'
#' @param model data.frame. Output of model-Functions.
#' @param data data.frame. Data input.
#' @param ts.ind tba
#' @param covariates tba
#' @param outcomes tba.
#' @return An object of class `data.frame`.
#' @export
#'
mlts_model_check <- function(model, data, ts.ind, covariates, outcomes){

  # a helper function to run a series of user input checks and print detailed warnings
  ## check if variable names are entered correctly
  check = !(ts.ind %in% colnames(data))
  if(sum(check) > 0){
    warning("Not all entries in ts.ind can be found in the data.")
  }
  check = !(names(covariates) %in% colnames(data))
  if(sum(check) > 0){
    warning("Not all entries in covariates can be found in the data.")
  }
  check = !(names(outcomes) %in% colnames(data))
  if(sum(check) > 0){
    warning("Not all entries in outcomes can be found in the data.")
  }

  ##



}