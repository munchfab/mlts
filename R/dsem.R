dsem <- function(formula = NULL, data,
                 y, id, beep,
                 iter = 2000, seed = NULL,
                 miss_handling = "remove") {

  dsem_terms <- terms(formula)
  mf <- model.frame(dsem_terms, data = data)

  return(dsem_terms)

  stan_data <- create_stan_data(
    data = data,
    y = y,
    id = id,
    beep = beep,
    miss_handling = miss_handling
  )

  # parameters to monitor
  pars <- c("btw_pred", # fixed effects of mu, ar, and (log) innovation variance
            "sigma", # random effect SDs
            "bcorr", # random effect correlations
            "bcov", # Var-Cov-matrix
            "y_rep")



  # draw samples
  # fit <- rstan::sampling(
  #   object = stanmodels$manifest_AR,
  #   chains = 4,
  #   cores = 4,
  #   iter = iter,
  #   data = stan_data,
  #   seed = seed,
  #   pars = pars
  # )

  # return(fit)

}
