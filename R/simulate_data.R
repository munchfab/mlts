sim_data_AR1_man <- function(N, TP,
                             mu, phi, logv,
                             sigma_mu, sigma_phi, sigma_logv,
                             mu_ec, sigma_ec,
                             cor_mu_phi, cor_mu_logv, cor_phi_logv,
                             cor_mu_ec, cor_phi_ec, cor_logv_ec,
                             miss = 0,
                             seed = NULL,
                             burnin = 100) {

  # create covariance matrix from inputs
  cov_mx <- diag(c(sigma_mu^2, sigma_phi^2, sigma_logv^2, sigma_ec^2))
  cov_mx[1, 2] <- cov_mx[2, 1] <- cor_mu_phi * sigma_mu * sigma_phi
  cov_mx[1, 3] <- cov_mx[3, 1] <- cor_mu_logv * sigma_mu * sigma_logv
  cov_mx[2, 3] <- cov_mx[3, 2] <- cor_phi_logv * sigma_phi * sigma_logv
  cov_mx[2, 4] <- cov_mx[4, 2] <- cor_phi_ec * sigma_ec * sigma_phi
  cov_mx[3, 4] <- cov_mx[4, 3] <- cor_logv_ec * sigma_ec * sigma_logv
  cov_mx[1, 4] <- cov_mx[4, 1] <- cor_mu_ec * sigma_ec * sigma_mu


  # create between-level data

  # Draw person parameters from multivariate normal distribution
  # To ensure stationarity, truncate phi between -1 and 1
  set.seed(seed)
  between <- data.frame(
    mnormt::rmtruncnorm(
      n = N,
      mean = c(mu, phi, logv, mu_ec),
      varcov = cov_mx,
      lower = c(-Inf, -1, -Inf, -Inf),
      upper = c(Inf, 1, Inf, Inf)
    )
  )
  names(between) <- c("mu", "phi", "logv", "ec")
  between$id <- rep(1:N)
  # return(between)

  # create within-level data in long format

  # create empty data frame
  within <- data.frame(
    id = rep(1:N, each = TP),
    y = NA
  )
  # generate within-level time series with individual AR1-effect
  for (i in 1:N) {
    ts <- as.numeric(
      arima.sim(
        # specify AR1-process with phi coefficient from between-level data
        model = list(ar = between[i, "phi"], order = c(1, 0, 0)),
        n = TP + burnin,
        rand.gen = rnorm, # innovations generated random normal
        sd = sqrt(exp(between[i, "logv"])) # with individual variance
      )
    )
    # add time series without burnin
    within[within$id == i, "y"] <- ts[(burnin + 1):length(ts)]
    # add time identifier
    within$TP <- rep(1:TP, times = N)
    # add trait level to time series
    within[within$id == i, "y"] <- within[within$id == i, "y"] + between$mu[i]
    # Introduce missings as stored in miss
    # N_miss <- between[between$id == i, "N_miss"]
    # N_complete <- TP - N_miss
    # within[within$id == i, "y"] = c(within[within$id == i, "y"][1:N_complete],
    #                                 rep(NA, times = N_miss))
  }
  return(within)
}
