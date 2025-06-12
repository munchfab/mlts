#' Generate Posterior Predictive Samples for Multilevel Latent Time Series Models
#'
#' @description
#' The `mlts_posterior_sample()` function generates replicated datasets from a fitted
#' \code{mlts} model using draws from the posterior distribution. The function can
#' simulate data under the population model or based on individual-specific (random effect) parameters.
#'
#' @param fit An object of class \code{mlts.fit}, as returned by a fitted model using \code{mlts()}.
#' @param draw_person_pars Logical. If \code{TRUE}, samples are generated using person-specific parameters (random effects).
#' If \code{FALSE}, only population-level parameters are used. Defaults to \code{FALSE}.
#' @param n_draws Integer. Number of posterior draws to use for simulating replicated datasets. Ignored if \code{draws} is provided. Defaults to 10.
#' @param draws Optional integer vector indicating specific posterior draw indices to use. If \code{NULL}, \code{n_draws} draws are randomly sampled from all available posterior samples.
#'
#' @details
#' The function extracts posterior samples of population-level (and optionally individual-level) parameters
#' from a fitted \code{mlts} model and simulates replicated datasets from the posterior predictive distribution.
#' Each replication corresponds to a different posterior draw and reflects uncertainty in the model's parameters.
#'
#' If \code{draw_person_pars = TRUE}, the function uses sampled person-specific random effects and covariate effects
#' from the posterior to generate new data at the individual level. This requires that the model was fitted with
#' \code{monitor_person_pars = TRUE} in \code{\link[mlts]{mlts_fit}}. If this condition is not met, the function will
#' throw an error.
#'
#' Posterior draws are either selected randomly (\code{n_draws}) or specified manually using the \code{draws} argument.
#' Optionally, left or right censoring is respected in the simulated data if such constraints were present in the model.
#'
#' @return A list of replicated datasets, each as a \code{data.frame} with columns:
#' \describe{
#'   \item{\code{Y_rep}}{Replication number.}
#'   \item{\code{ID}}{Subject/cluster ID.}
#'   \item{\code{time}}{Time point.}
#'   \item{...}{One column per time-series variable defined in the model.}
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate 20 replications from the posterior
#' y_reps <- mlts_posterior_sample(fit = my_model_fit, n_draws = 20)
#'
#' # Include person-specific parameters in simulation
#' y_reps <- mlts_posterior_sample(fit = my_model_fit, draw_person_pars = TRUE)
#'
#' # Use specific posterior draws
#' y_reps <- mlts_posterior_sample(fit = my_model_fit, draws = c(10, 50, 100))
#' }
#'
#' @seealso \code{\link[mlts]{mlts_pp_check}} for plotting posterior predictive checks.
#' @export

mlts_posterior_sample <- function(
    fit,
    draw_person_pars = FALSE,
    n_draws = 10,
    draws = NULL
){

  # get model infos
  infos = mlts_model_eval(fit$model)

  # get number of posterior draws
  n_mcmc = prod(dim(fit$posteriors[,,1]))

  # obtain posterior draws of model parameters
  mm_pars = fit$pop.pars.summary[fit$pop.pars.summary$Model == "Measurement",]
  if(infos$isLatent == TRUE){
    mm_samples = rstan::extract(fit$stanfit, pars = mm_pars$Param_stan)
  }
  cor_pars = fit$pop.pars.summary[fit$pop.pars.summary$Param_Label == "Innovation correlation",]
  if(infos$n_inno_cors > 1){
    cor_samples = rstan::extract(fit$stanfit, pars = cor_pars$Param_stan)
  }

  if(draw_person_pars == TRUE){
    W = fit$standata$W
    gamma_pars = fit$pop.pars.summary[grepl(fit$pop.pars.summary$Param_stan, pattern = "gammas"),]
    gamma_samples = rstan::extract(fit$stanfit, pars = "gammas")
    sd_R_samples = rstan::extract(fit$stanfit, pars = "sd_R")
    bcorr_samples = rstan::extract(fit$stanfit, pars = "bcorr")
    if(infos$n_cov > 1){
      re_pred_pars = fit$pop.pars.summary[grepl(fit$pop.pars.summary$Param_stan, pattern = "b_re_pred"),]
      re_pred_samples = rstan::extract(fit$stanfit, pars = re_pred_pars$Param_stan)
    }
  }

  # select draws
  if(is.null(draws)){draws_use = sample(x = 1:n_mcmc, size = n_draws)}

  # list of replications
  y_reps = list()

  # check if posterior samples of person parameters exist
  if(draw_person_pars == FALSE & is.na(fit$person.pars.summary)[1]){
    stop("Posterior samples of person-specific parameters are not available.
         Consider setting `monitor_person_pars = TRUE` in `mlts_fit`.")
  }

  for(i in 1:n_draws){

    # obtain posterior values of specific draw
    if(infos$isLatent == TRUE){
      mm_pars$sample = NA
      for(j in 1:nrow(mm_pars)){
        mm_pars$sample[j] = mm_samples[[mm_pars$Param_stan[j]]][draws_use[i]]
      }
    }

    if(infos$n_inno_cors > 1){
      cor_pars$sample = NA
      for(j in 1:nrow(cor_pars)){
        cor_pars$sample[j] = cor_samples[[cor_pars$Param_stan[j]]][draws_use[i]]
      }
    }

    # sample cluster under the population model
    if(draw_person_pars == TRUE){
      btw = get_new_person_par_mat(
        fit = fit, infos = infos, iter = draws_use[i], gamma_pars = gamma_pars,
        gamma_samples = gamma_samples, sd_R_samples = sd_R_samples,
        bcorr_samples = bcorr_samples, re_pred_pars = re_pred_pars,
        re_pred_samples = re_pred_samples, W = W
      )
    }

    # use posterior samples
    if(draw_person_pars == FALSE){
      btw = get_person_par_mat(fit, infos = infos, i)
    }

    # use observed exogenous time series variables
    exogenous = NULL
    if (any(infos$is_wcen==0)) {
      if( infos$isLatent == FALSE ){
        exogenous = as.matrix(fit$data[,fit$standata$ts[infos$is_wcen==0]])
      } else {
        exogenous = as.matrix(fit$data[,fit$standata$ts[infos$p_is_wcen==0]])
      }
    }

    # sample time series data
    reps = mlts_sim_within(
      infos = infos,
      burn.in = 0,
      N = fit$standata$N,
      TP = fit$standata$N_obs_id,
      btw = btw,
      mm_pars = mm_pars,
      cor_pars = cor_pars,
      exogenous = exogenous)

    # add proper names
    colnames(reps) <- c("ID", "time", fit$standata$ts)

    # add censoring
    if(!is.null(attr(fit$model, which = "censor_left"))){
      reps[,fit$standata$ts][reps[,fit$standata$ts] < fit$standata$censL_val] <- fit$standata$censL_val
    }

    if(!is.null(attr(fit$model, which = "censor_right"))){
      reps[,fit$standata$ts][reps[,fit$standata$ts] > fit$standata$censR_val] <- fit$standata$censR_val
    }

    # add replication number for plotting
    y_reps[[i]] <- cbind(
      "Y_rep" = i,
      reps
    )

  }

  return(y_reps)

}
