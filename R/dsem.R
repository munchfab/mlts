#' Title
#'
#' @param formula
#' @param data
#' @param y
#' @param id
#' @param beep
#' @param iter
#' @param seed
#' @param miss_handling
#'
#' @return
#' @export
#'
#' @examples
dsem <- function(formula, data,
                 y = NULL, id = NULL, beep = NULL,
                 iter = 2000, seed = NULL,
                 miss_handling = "remove") {

  # parse formula for variables to be passed to Stan
  dsem_terms <- dsem_formula(formula)

  # create Stan data from parsed formula
  stan_data <- create_stan_data(
    data = data,
    y = dsem_terms$response,
    id = dsem_terms$id,
    beep = dsem_terms$beep,
    miss_handling = miss_handling
  )

  # parameters to monitor
  # needs updating for different models
  pars <- c("btw_pred", # fixed effects of mu, ar, and (log) innovation variance
            "sigma", # random effect SDs
            "bcorr", # random effect correlations
            "bcov", # Var-Cov-matrix
            "y_rep")


  # draw samples for respective model
  # if latent variables are present, use latent variable model
  if (grepl("lv", specials)) {
    NULL
  } else if (grepl("ar", specials)) {
    fit <- rstan::sampling(
      object = stanmodels$manifest_AR,
      chains = 4,
      cores = 4,
      iter = iter,
      data = stan_data,
      seed = seed,
      pars = pars
    )
  }

  return(fit)
}
