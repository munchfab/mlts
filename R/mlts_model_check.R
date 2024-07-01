#' Helper function to run a series of user input checks and print detailed warnings
#'
#' @param model `data.frame`. Output of model-Functions.
#' @param data `data.frame`. Data input.
#' @param ts tba
#' @param covariates tba
#' @param outcomes tba.
#' @return An object of class `data.frame`.
#' @noRd
#'
mlts_model_check <- function(model, data, ts, covariates, outcomes){

  # check if variable names are entered correctly
  check = !(ts %in% colnames(data))
  if(sum(check) > 0){
    warning("Not all entries in ts can be found in the data.")
  }
  check = !(names(covariates) %in% colnames(data))
  if(sum(check) > 0){
    warning("Not all entries in covariates can be found in the data.")
  }
  check = !(names(outcomes) %in% colnames(data))
  if(sum(check) > 0){
    warning("Not all entries in outcomes can be found in the data.")
  }
}
