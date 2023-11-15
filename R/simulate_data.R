sim_data_AR1_man <- function(N, TP,
                             mu, phi, logv,
                             sigma_mu, sigma_phi, sigma_logv,
                             mu_ec, sigma_ec,
                             cor_mu_phi, cor_mu_logv, cor_phi_logv,
                             cor_mu_ec, cor_phi_ec, cor_logv_ec,
                             miss = 0,
                             burnin = 100) {
  # create covariance matrix from inputs
  cov_mx <- diag(c(sigma_mu^2, sigma_phi^2, sigma_logv^2, sigma_ec^2))
  cov_mx[1, 2] <- cov_mx[2, 1] <- cor_mu_phi * sigma_mu * sigma_phi
  cov_mx[1, 3] <- cov_mx[3, 1] <- cor_mu_logv * sigma_mu * sigma_logv
  cov_mx[2, 3] <- cov_mx[3, 2] <- cor_phi_logv * sigma_phi * sigma_logv
  cov_mx[2, 4] <- cov_mx[4, 2] <- cor_phi_ec * sigma_ec * sigma_phi
  cov_mx[3, 4] <- cov_mx[4, 3] <- cor_logv_ec * sigma_ec * sigma_logv
  cov_mx[1, 4] <- cov_mx[4, 1] <- cor_mu_ec * sigma_ec * sigma_mu


  # create between-level-data #################################################
  # Draw person parameter from multivariate normal distribution
  # To ensure stationarity, truncate phi between -1 and 1
  btw <- data.frame(
    mnormt::rmtruncnorm(
      n = N,
      mean = c(mu, phi, logv, mu_ec),
      varcov = cov_mx,
      lower = c(-Inf, -1, -Inf, -Inf),
      upper = c(Inf, 1, Inf, Inf)
    )
  )

}
