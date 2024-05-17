#' Get Standardized Estimates for an mlts Model
#'
#' @param object `mltsfit`. Output of \code{\link[mlts]{mlts_model}} and related functions.
#' @param what character. Get between-level standardized estimates (`what = "between"`, the default),
#' within-level standardized estimates averaged over clusters (`what = "within"`), or both (`what = "both"`).
#' @param prob A value between 0 and 1 to indicate the width of the credible
#' interval. Default is .95.
#' @param digits Number of digits. Default is 3.
#' @param add_cluster_std logical. If `what = "within"`, within-level standardized effects for each cluster
#' are included in the output (defaults to `FALSE`).
#' @return A `list` containing between- and within-level standardized parameters.
#' @export
#'
#' @examples
#' \donttest{
#' # build simple vector-autoregressive mlts model for two time-series variables
#' var_model <- mlts_model(q = 2)
#'
#' # fit model with (artificial) dataset ts_data
#' fit <- mlts_fit(
#'   model = var_model,
#'   data = ts_data,
#'   ts = c("Y1", "Y2"), # time-series variables
#'   id = "ID", # identifier variable
#'   time = "time", # time variable
#'   tinterval = 1, # interval for approximation of continuous-time dynamic model,
#'   monitor_person_pars = TRUE # person parameters need to be sampled for standardization
#' )
#'
#' # inspect standardized parameter estimates
#' mlts_standardized(fit)
#' }
#'
mlts_standardized <- function(object, what = c("between", "within", "both"),
                              digits = 3, prob = .95, add_cluster_std = FALSE
){

  what <- match.arg(what)
  result <- list()

  # get model infos
  infos <- mlts_model_eval(object$model)

  # make sure object is of class mltsfit
  # if(class(object) != "mltsfit")
  if (!inherits(object, "mltsfit")) {
    stop("Input of `object` should be of class 'mltsfit'.")
  }

  # run between-level standardization
  if(what == "between" | what == "both"){
  std_btw = mlts_standardized_btw(object = object, digits = digits, prob = prob)
  result[["Between-level standardized"]] = std_btw
  }


  if(what == "within" | what == "both"){
  if(length(object$person.pars.summary) == 1){
    stop("\n To obtain standardized estimates per cluster and the average standardized estimate(s) across clusters, refit the model with `monitor_person_pars = TRUE`.",
         "\n Note that for multiple indicator models (p > 1), `get_SD_latent` needs to be set to TRUE.")
  }

  if(infos$isLatent == TRUE & object$standata$standardized == 0){
    stop("\n Variance of latent factor scores not available for standardization
            of dynamic model parameters. Refit the model with get_SD_latent = TRUE
            to obtain standardized estimates.")
  }

    # run within-level standardization
    std_within = mlts_standardize_within(object = object, digits = digits, prob = prob,
                            add_cluster_std = add_cluster_std)
    result[["Within-level standardized effects averaged over clusters"]] =
      std_within$`Within-level standardidzed effects averaged over clusters`
    if(add_cluster_std == TRUE){
      result[["Within-level standardized effects by cluster"]] =
        std_within$`Within-level standardized effects by cluster`
    }
  }

  return(result)

}
