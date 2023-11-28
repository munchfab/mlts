# Fit stan model using artificial AR1 data set

devtools::load_all()

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


# fit model
fit <- dsem(
  formula = y ~ TP + (TP | id),
  data = ar1_data,
  y = "y", id = "id", beep = "TP",
  iter = 50, seed = 12
)

print(fit)

# ar1_stan_data <- create_stan_data(
#   data = ar1_data, id = "id", beep = "TP", y = "y",
#   pred_random = NULL,
#   outcome = NULL,
#   out_predictors = NULL,
#   standardize_out = TRUE,
#   overnight_lags = NULL,
#   miss_handling = "remove",
#   random_innovations = TRUE,
#   add_mplus_data = TRUE
# )

# load artifical data set
# load("./data/ar1_stan_data.rda")

# compile model
# not needed if precompiled
# ar1_model <- rstan::stan_model(
#   file = "./src/stan_files/manifest_AR.stan",
#   model_name = "manifest_AR"
# )

# # parameters to monitor
# pars <- c("btw_pred", # fixed effects of mu, ar, and (log) innovation variance
#           "sigma", # random effect SDs
#           "bcorr", # random effect correlations
#           "bcov", # Var-Cov-matrix
#           "y_rep")
#
# # draw samples
# ar1_fit <- rstan::sampling(
#   object = stanmodels$manifest_AR,
#   chains = 4,
#   cores = 4,
#   iter = 3000,
#   data = ar1_stan_data,
#   seed = 1234,
#   pars = pars
# )
#
# # print summary
# print(ar1_fit, pars = pars[!pars == "y_rep"])
# print(ar1_fit, pars = "y_rep")

# extract posterior predictions
y_rep <- as.matrix(
  ar1_fit,
  pars = paste0("y_rep[", 1:ar1_stan_data$N_obs, "]"))[sample(500), ]


# posterior predictive checking: density
ppc_dens <- bayesplot::ppc_dens_overlay(
  ar1_stan_data$y,
  y_rep
)
# view
ppc_dens
