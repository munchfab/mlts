# manifest
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

# latent
sim_data_AR1_lat <- function(N, TP, IND,
                             lambda_w, sigma_w,
                             lambda_b, alpha, sigma_b,
                             mu, phi, logv,
                             sigma_mu, sigma_phi, sigma_logv,
                             mu_ec, sigma_ec,
                             cor_mu_phi, cor_mu_logv, cor_phi_logv,
                             cor_mu_ec, cor_phi_ec, cor_logv_ec,
                             miss = 0,
                             seed = NULL,
                             burnin = 100) {

  # mus
  bmus <- c(mu, phi, logv, mu_ec)

  # create covariance matrix from inputs
  cov_mx <- diag(c(sigma_mu^2, sigma_phi^2, sigma_logv^2, sigma_ec^2))

  # within-level
  lambda_w <- lambda_w # factor loadings
  sigma_w <- diag(sigma_w) # covariance matrix of measurement errors

  # create between-level data
  lambda_b <- lambda_b
  alphas <- matrix(alpha, ncol = 1) # item intercepts
  sigma_b <- diag(sigma_b)   # covariance matrix of measurement errors

  # draw person parameter from multivariate normal distribution
  btw <- mnormt::rmnorm(n = N, mean = bmus, varcov = cov_mx)

  # generate latent process
  eta_w <- list()

  for (pp in 1:N) {
    eta_w[[pp]] <- c(1:TP)
    eta_w[[pp]][1] <- rnorm(1, 0, sqrt(exp(btw[pp,3])))
    for (tt in 2:TP) {
      eta_w[[pp]][tt] <- btw[pp, 2] * eta_w[[pp]][tt - 1] + rnorm(1, 0, sqrt(exp(btw[pp, 3])))
    }
  }

  # check if data are reasonable
  # psych::autoR(etaW[[1]])
  # btw[1,2]
  # summary(lm(eta_w[[1]][2:TP] ~ eta_w[[1]][1:(TP - 1)]))$sigma
  # sqrt(exp(btw[1, 3]))

  # generate data
  y <- array(NA, dim = c(N, TP, IND))

  for(pp in 1:N) {
    # sample error vectors
    # within
    eps_w <- mnormt::rmnorm(n = TP, rep(0, IND), sigma_w) # sample error vector
    y_w <- matrix(NA, ncol = IND, nrow = TP)
    y_w <- matrix(eta_w[[pp]], ncol = 1) %*% t(cbind(lambda_w)) + eps_w

    eps_b <- mnormt::rmnorm(1, rep(0, IND), sigma_b) # sample error vector
    y_b <- matrix(NA, ncol = IND, nrow = 1)
    y_b <- t(alphas) + matrix(btw[pp, 1], ncol = 1) %*% t(cbind(lambda_b)) + eps_b

    y[pp, , ] <- matrix(y_b, nrow = TP, ncol = IND, byrow = T) + y_w
  }

  # retransform data
  y_new <- array(data = 0, dim = c(IND, N, TP))

  for (ii in 1:IND) {
    for (pp in 1:N) {
      for (tt in 1:TP) {
        y_new[ii, pp, tt] <- y[pp, tt, ii]
      }
    }
  }

  # bring data into long format data frame
  df <- data.frame(
    id = rep(1:N, each = TP),
    TP = rep(1:TP, times = N)
  )

  for (ii in 1:IND) {
    for (pp in 1:N) {
      for (tt in 1:TP) {
        df[df$id == pp & df$TP == tt, 2 + ii] <- y_new[ii, pp, tt]
      }
    }
  }

  # check if retranformation was correct
  # mean(df$V3[df$N == 2])
  # mean(y_new[1, 2, ])
  # mean(df$V5[df$TP == 5])
  # mean(y_new[3, , 5])

  # proper column names
  colnames(df) <- c("id", "TP", paste0("y", 1:IND))

  # return data and parameter inputs
  return(df)
}

