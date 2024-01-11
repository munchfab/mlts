#' Title
#'
#' @param y character. The variable in `data` that contains the observations
#' of the time-varying construct.
#' @param id character. The variable in `data` that identifies the person or
#' observational unit.
#' @param beep character. The variable in `data` that contains the running
#' beep from 1 to TP for each person.
#' @param data An object of class `data.frame` (or one that can be coerced
#' to that class) in long format.
#' @param miss_handling character. Should missings be removed (`"remove"`) or
#' imputed (`"impute"`)?
#' @param iter Number of total iterations per chain (including warmup; defaults
#' to 2000).
#' @param seed The seed for random number generation to make results
#' reproducible. If \code{NULL} (the default), \pkg{Stan} will set the seed
#' randomly.
#'
#' @return An object of class `stanfit`.
#' @export
#'
#' @examples
dsem_ar <- function(y, id, beep, data,
                    miss_handling = "remove",
                    iter = 2000, seed = NULL) {

  # VarModelBuild

  # create Stan data list
  stan_data <- create_stan_data(
    data = data,
    y = y,
    id = id,
    beep = beep,
    miss_handling = miss_handling
  )

  # parameters to monitor
  # needs updating for different models
  pars <- c("btw_pred", # fixed effects of mu, ar, and (log) innovation variance
            "sigma", # random effect SDs
            "bcorr", # random effect correlations
            "bcov", # Var-Cov-matrix
            "y_rep")

  # draw samples from model
  fit <- rstan::sampling(
    object = stanmodels$manifest_AR,
    chains = 4,
    cores = 4,
    iter = iter,
    data = stan_data,
    seed = seed,
    pars = pars
  )

  return(fit)
}
