#' Get Standardized Estimates for Between-level Regression Parameters of an mlts Model
#'
#' @param object `mltsfit`. Output of \code{\link[mlts]{mlts_model}} and related functions.
#' @param prob A value between 0 and 1 to indicate the width of the credible
#' interval. Default is .95.
#' @param digits Number of digits.
#' @return An object of class `list`.
#' @export
#'
#' @examples
#' \donttest{
#' # build simple vector-autoregressive mlts model for two time-series variables
#' var_model <- mlts_model(
#'     q = 2,
#'     ranef_pred = c("covariate1", "covariate2"),
#'     out_pred = c("outcome1", "outcome2")
#' )
#'
#' # fit model with (artificial) dataset ts_data
#' fit <- mlts_fit(
#'   model = var_model,
#'   data = ts_data,
#'   ts = c("Y1", "Y2"), # time-series variables
#'   id = "ID", # identifier variable
#'   tinterval = 1 # interval for approximation of continuous-time dynamic model,
#' )
#'
#' # inspect standardized parameter estimates
#' mlts_standardized_btw(fit)
#' }
#'
mlts_standardized_btw <- function(object, digit = 3, prob = .95
){

  # make sure object is of class mltsfit
  if(class(object) != "mltsfit"){
    stop("Input of `object` should be of class 'mltsfit'.")
  }

  # get model infos
  infos <- mlts_model_eval(object$model)

  # get CIs as specified by user
  alpha = 1 - prob
  probs = c(alpha/2, 1-alpha/2)
  prob.cols = paste0(c(100-(100-(alpha/2)*100), 100-(alpha/2)*100), "%")
  # result columns
  result.cols = c("Std.Est", "SD", prob.cols)
  result.list <- list() # final object
  result <- data.frame()

  # Between-level standardized effects =========================================
  ## Get posteriors of all random effects by parameter
  chains <- as.numeric(object$stanfit@sim$chains)
  warmup <- as.numeric(object$stanfit@sim$warmup)
  iter <- as.numeric(object$stanfit@sim$iter) -warmup
  n_random <- object$standata$n_random
  re_sds = array(dim = c(chains, iter, n_random))
  # get SDs of RE predictor variables
  Wvars_sds = apply(object$standata$W, MARGIN = 2, FUN = sd)

  # get/calculate model-implied RE variances ----
  ## get labels of (residual) random effect SDs
  labs_re = paste0("sd_R[",1:n_random,"]")
  sigma_RE = rstan::extract(object$stanfit, pars = labs_re)
  # number of covariates used for prediction of random effects
  n_cov = object$standata$n_cov
  # as default use sigma of random effects and update respectively
  Var_RE = sigma_RE
  for(i in 1:n_random){
    Var_RE[[i]] = Var_RE[[i]]^2
  }

  # Part 1: Standardized estimates for random effects regressed on covariate(s)
  if(object$standata$n_cov > 1){
    for(i in 1:object$standata$n_cov_bs){
      re_pos.i = object$standata$n_cov_mat[i,2]
      # get unstandardized estimates of re prediction paramerters
      b <- rstan::extract(object$stanfit, pars = paste0("b_re_pred[",i,"]"))
      # update RE variances based on prediction parameter contributions
      Var_RE[[re_pos.i]] = Var_RE[[re_pos.i]] + b[[1]]^2 * Wvars_sds[object$standata$n_cov_mat[i,1]]^2
    }

    # now use SDs of RE to standardize regression parameters:
    # prepare object to store results
    re_pred_std = infos$RE.PREDS[, c("Type", "Param")]
    re_pred_std[,result.cols] = NA
    for(k in 1:nrow(infos$RE.PREDS)){

      # get SD of covariate
      sd_x = Wvars_sds[object$standata$n_cov_mat[k,1]]
      # get position of RE on RE_par_SD
      sd_y = sqrt(Var_RE[[object$standata$n_cov_mat[k,2]]])
      # get parameter label in stan model
      b <- rstan::extract(object$stanfit, pars = paste0("b_re_pred[",k,"]"))
      b_std <- b[[1]] * sd_x / sd_y
      # save summary statistics
      re_pred_std[k, result.cols] = round(c(
        mean(unlist(b_std)),
        sd(unlist(b_std)),
        quantile(unlist(b_std), c(probs))),digits = digit)
    }

    result = rbind(result, re_pred_std)
  }

  # Part 2: Standardized estimates for outcomes regressed on random effects
  if(object$standata$n_out > 0){
    # prepare final object to store results
    out_pred_std = infos$OUT[, c("Type", "Param")]
    out_pred_std[,result.cols] = NA

    # get SDs of outcomes
    SDs_out = apply(object$standata$out, FUN = sd, MARGIN = 1)
    # add variances of additional covariates used as predictor
    if(object$standata$n_z > 0){
      for(i in 1:object$standata$n_z){
        Var_RE[[n_random+i]] <- var(object$standata$Z[,i])
      }
    }

    for(i in 1:nrow(infos$OUT)){
      # get unstandardized estimate
      lab_out = object$pop.pars.summary$Param_stan[which(object$pop.pars.summary$Param == infos$OUT$Param[i])]
      b = rstan::extract(object$stanfit, pars = c(lab_out))
      sd_y = SDs_out[infos$OUT$out_var_no[i]]
      sd_x = sqrt(Var_RE[[infos$OUT$Pred_no[i]]])
      b_std = b[[1]] * sd_x / sd_y
      # save summary statistics
      out_pred_std[i, result.cols] = round(c(
        mean(unlist(b_std)),
        sd(unlist(b_std)),
        quantile(unlist(b_std), c(probs))),digits = digit)
    }

    row.names(result) <- NULL
    result = rbind(result, out_pred_std)
  }

  result.list[["Between-level standardized"]] <- result


  return(result.list)
}
