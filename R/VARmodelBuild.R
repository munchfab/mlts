#' Title
#'
#' @param q integer. The number of time-varying constructs.
#' @param p integer. For multiple-indicator models, specify a vector of length
#' `q` with the number of manifest indicators per construct.
#' @return An object of class `data.frame`.
#' @export
#'
VARmodelBuild <- function(q, p = NULL){

  # checks needed?
  ### Müssen wir hier einen manuellen Check einbauen, ob p entweder
  ### length(p) == 1 oder length(p) == q ist?
  # setze ich in testthat um
  if(length(p) == 1){
    p = rep(p, times = q)
    if (q > 1) {
      message("Note: The number of indicators is assumed to be ", p[1],
              " for each latent variable. If this is not intended, please",
              " specify a vector of length q containing the number of",
              " indicators for each latent construct",
              " (see Vignettes for examples).")
    }
  }
  if (length(p) != q) {
    stop("If multiple indicators are used for latent constructs,",
         " p should be a vector of length q containing the number of",
         " indicators for each latent construct (see Vignettes for examples).")
  }


  # Structural Model ===========================================================
  n_mus = q                                 # trait level parameters
  mus_pars = paste0("mu_",1:n_mus)

  n_phi = q^2                               # dynamic parameters
  phi_pars = paste0("phi_",rep(1:q, each = q), rep(1:q, times = q))

  n_sigma = q                               # innovation variances
  sigma_pars = paste0("ln.sigma2_", 1:q)

  n_covs = (q *(q-1)) / 2                   # innovation covariances
  qs = c()
  ps = c()
  for(i in 1:(q-1)){
    qs = c(qs, rep(i, each = q-i))
  }
  for(i in 2:(q)){
    ps = c(ps, rep(i:q, 1))
  }
  cov_pars = paste0("ln.sigma_", qs, ps)

  # ---

  if(q > 1){

message("Note: When specifying a VAR(1) model with person-specific innovation
covariances, a latent-variable approach will be used which affords introducing
contraints on the loading parameters of the latent covariance factor(s).
(see Hamaker et al., 2018). If innovation covariances are set as a constant
or set to 0 in a subsequent step, this warning can be ignored.")

    pars = c(mus_pars, phi_pars, sigma_pars, cov_pars)

    # combine a information in a data frame
    FE = data.frame(
      "Model" = "Structural",
      "Level" = "Within",
      "Type" = "Fix effect",
      "Param" = pars,
      "Param_Label" = c(
        rep("Trait", n_mus),
        rep("Dynamic", n_phi),
        rep("Log Innovation Variance", n_sigma),
        rep("Log Innovation Covariance", n_covs)
      ),
      "isRandom" = 1
    )
    RE = data.frame(
      "Model" = "Structural",
      "Level" = "Between",
      "Type" = "Random effect SD",
      "Param" = paste0("sigma_",pars),
      "Param_Label" = c(
        rep("Trait", n_mus),
        rep("Dynamic", n_phi),
        rep("Log Innovation Variance", n_sigma),
        rep("Log Innovation Covariance", n_covs)
      ),
      "isRandom" = 0
    )
    ## combine
    df.pars = rbind(FE, RE)
  } else {
    pars = c(mus_pars, phi_pars, sigma_pars)

    # combine a information in a data frame
    ## fixed effects
    FE = data.frame(
      "Model" = "Structural",
      "Level" = "Within",
      "Type" = "Fix effect",
      "Param" = pars,
      "Param_Label" = c(
        rep("Trait", n_mus),
        rep("Dynamic", n_phi),
        rep("Log Innovation Variance", n_sigma)
      ), "isRandom" = 1
    )
    ## random effect variances
    RE = data.frame(
      "Model" = "Structural",
      "Level" = "Between",
      "Type" = "Random effect SD",
      "Param" = paste0("sigma_",pars),
      "Param_Label" = c(
        rep("Trait", n_mus),
        rep("Dynamic", n_phi),
        rep("Log Innovation Variance", n_sigma)
      ),
      "isRandom" = 0
    )
  }

    ## add random effect correlations
    rand.pars = FE[FE$isRandom == 1,"Param"]
    n_rand = length(rand.pars)
    btw.cov_pars = c()
    if(n_rand>1){
      n_cors = (n_rand * (n_rand-1))/2
      qs = c()
      ps = c()
      for(i in 1:(n_rand-1)){
        qs = c(qs, rep(rand.pars[i], each = n_rand-i))
      }
      for(i in 2:n_rand){
        ps = c(ps, rep(rand.pars[i:n_rand], 1))
      }

      btw.cov_pars = paste0("r_", qs,".", ps)

      ## random effect correlations
      REcors = data.frame(
        "Model" = "Structural",
        "Level" = "Between",
        "Type" = rep("RE correlation",n_cors),
        "Param" = btw.cov_pars,
        "Param_Label" = "RE Cor",
        "isRandom" = 0
      )


    ## combine
    VARmodel = rbind(FE, RE, REcors)
  }

  # ADD DEFAULT PRIORS ========================================================
  VARmodel = VARmodelPriors(VARmodel = VARmodel, default = T)

  # MEASUREMENT MODEL =========================================================
  if(!is.null(p)){
    VARmodel = VARmodelMeasurement(VARmodel = VARmodel, q = q, p = p)
  }

  return(VARmodel)

}
