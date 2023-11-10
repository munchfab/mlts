# Fit stan model using artificial AR1 data set

# load artifical data set
load("./data/ar1_stan_data.rda")

# compile model
ar1_model <- rstan::stan_model(
  file = "./src/stan_files/manifest_AR.stan",
  model_name = "manifest_AR"
)

# parameters to monitor
pars <- c("btw_pred", # fixed effects of mu, ar, and (log) innovation variance
          "sigma", # random effect SDs
          "bcorr", # random effect correlations
          "bcov") # Var-Cov-matrix

# draw samples
ar1_fit <- rstan::sampling(
  object = ar1_model,
  chains = 4,
  cores = 4,
  iter = 3000,
  data = ar1_stan_data,
  seed = 1234,
  pars = pars
)

# print summary
print(ar1_fit)
