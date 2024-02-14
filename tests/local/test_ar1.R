# Fit stan model using artificial AR1 data set

## manifest model =============================================================

devtools::load_all()

ar1_data <- VARmodelSim(VARmodel = VARmodel)

# create artificial data
ar1_data <- sim_data_AR1_man(
  N = 100,
  TP = 70,
  mu = 10, phi = .3, logv = 2,
  sigma_mu = 1.4, sigma_phi = .2, sigma_logv = .7,
  mu_ec = 5, sigma_ec = 2,
  cor_mu_phi = .3, cor_mu_logv = .3, cor_phi_logv = .3,
  cor_mu_ec = .3, cor_phi_ec = .3, cor_logv_ec = .3,
  seed = 1234
)

ar1_data_stan <- create_stan_data(
  data = ar1_data,
  y = "y", id = "id", beep = "TP",
  miss_handling = "remove"
)


# fit model
ar1_fit <- dsem_ar(
  y = "y", id = "id", beep = "TP",
  data = ar1_data,
  # y = "y", id = "id", beep = "TP",
  iter = 2000, seed = 12
)

# summary
print(ar1_fit)

# extract posterior predictions
y_rep <- as.matrix(
  ar1_fit,
  pars = paste0("y_rep[", 1:length(ar1_data$id), "]"))[sample(500), ]

# posterior predictive checking: density
ppc_dens <- bayesplot::ppc_dens_overlay(
  ar1_data$y,
  y_rep
)

# view
ppc_dens



## latent model ===============================================================

devtools::load_all()

# create artificial data
ar1_data_lat <- sim_data_AR1_lat(
  N = 100,
  TP = 70,
  IND = 3,
  lambda_w = c(.8, .8, .7), sigma_w = c(.3, .2, .3),
  lambda_b = c(.8, .8, .7), alpha = c(2, 2, 2),
  sigma_b = c(.4, .3, .2),
  mu = 10, phi = .3, logv = 2,
  sigma_mu = 1.4, sigma_phi = .2, sigma_logv = .7,
  mu_ec = 5, sigma_ec = 2,
  cor_mu_phi = .3, cor_mu_logv = .3, cor_phi_logv = .3,
  cor_mu_ec = .3, cor_phi_ec = .3, cor_logv_ec = .3,
  seed = 1234
)

devtools::load_all()



# create stan data
stan_data <- create_stan_data(
  data = ar1_data_lat,
  id = "id", beep = "TP",
  y = c("y1", "y2", "y3"),
  # N_ind = 3,
  pred_random = NULL,
  outcome = NULL,
  out_predictors = NULL,
  standardize_out = TRUE,
  overnight_lags = NULL,
  miss_handling = "remove",
  random_innovations = TRUE,
  add_mplus_data = TRUE
)


# load artifical data set
# load("./data/ar1_stan_data.rda")

# compile model
# not needed if precompiled
ar1_model <- rstan::stan_model(
  file = "./inst/stan/latent_AR.stan",
  model_name = "latent_AR"
)


# parameters to monitor
pars <- c("btw_pred", # fixed effects of mu, ar, and (log) innovation variance
          "sigma", # random effect SDs
          "bcorr", # random effect correlations
          "bcov", # Var-Cov-matrix
          "yB_rep"
          # "yW_rep"
          )
#
# draw samples
ar1_fit <- rstan::sampling(
  object = ar1_model,
  chains = 4,
  cores = 4,
  iter = 3000,
  data = stan_data,
  seed = 1234,
  pars = pars
)

print(ar1_fit)

rstan::traceplot(ar1_fit, pars = pars)
#
# # print summary
# print(ar1_fit, pars = pars[!pars == "y_rep"])
# print(ar1_fit, pars = "y_rep")

