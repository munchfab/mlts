#' Title
#'
#' @param q integer. The number of time-varying constructs.
#' @param p integer. For multiple-indicator models, specify a vector of length
#' `q` with the number of manifest indicators per construct.
#' @param FiVARmodeledInnoVars logical. FiVARmodel all random effect variances (eVARmodelcept those
#' of individual traits) to zero.
#' @param FiVARmodeledCovs logical. Set all innovation covariances to a constant value.
#' @param CovsZero logical. Set to TRUE to treat all innovations as independent.
#' @param FEis0 character. A character vector to indeVARmodel which fiVARmodeled model parameters
#' should be fiVARmodeled to zero (Note: this results in removing the random effect
#' variance of the respective parameter).
#' @param REis0 logical. Set to TRUE to treat all innovations as independent.
#' @param btw.factor Logical. If `TRUE` (the default), a common between-level factor
#' is modeled across all indicator variables. If `FALSE`, instead of a between-level
#' factor, indicator mean levels will be included as individual (random) effects stemming
#' from a joint multivariate normal distribution.
#' @param btw.model A list to indicate for which manifest indicator variables a common
#' between-level factor should be modeled (see Details for detailed instructions).
#' At this point restricted to one factor per latent construct.
#' @param RE.pred character vector or a named list. Include between-level covariate(s)
#' as predictor(s) of all random effects in `VARmodel` by entering a vector of unique variable
#' names. Alternatively, to include between-level covariates or differing sets of
#' between-level covariates as predictors of specific random effects, a named
#' list (using the `param`-labels in `VARmodel`) can be entered (see details).
#' Note that if a named list is provided, all names that do not match random
#' parameters in `VARmodel` will be ignored.
#' @param out.pred character vector or a named list. Include between-level covariate(s)
#' as predictor(s) of all random effects in `VARmodel` by entering a vector of unique variable
#' names. Alternatively, to include between-level covariates or differing sets of
#' between-level covariates as predictors of specific random effects, a named
#' list (using the `param`-labels in `VARmodel`) can be entered (see details).
#' @return An object of class `data.frame`.
#' @export
#'
mlts_model <- function(q, p = NULL, maxLag = c(1,2,3),
                          btw.factor = TRUE, btw.model = NULL,
                          FixDynamics = F, FixInnoVars = F,
                          FixInnoCovs = F, InnoCovsZero = NULL,
                          FEis0 = NULL, REis0 = NULL,
                          RE.pred = NULL, out.pred=NULL, out.pred.add.btw = NULL){

  if(length(maxLag) == 3){
    maxLag = 1
  }

  if(length(p) == 1){
    p = rep(p, times = q)
    if (q > 1) {
      warning("Note: The number of indicators is assumed to be ", p[1],
              " for each latent variable. If this is not intended, please",
              " specify a vector of length q containing the number of",
              " indicators for each latent construct",
              " (see Vignettes for examples).")
    }
  }
  if (length(p) != q & !is.null(p)) {
    stop("If multiple indicators are used for latent constructs,",
         " p should be a vector of length q containing the number of",
         " indicators for each latent construct (see Vignettes for examples).")
  }


  # Structural Model ===========================================================
  n_mus = q                                 # trait level parameters
  mus_pars = paste0("mu_",1:n_mus)

  n_phi = (q^2)*maxLag                        # dynamic parameters
  phi_order = rep(paste0("phi(",1:maxLag,")_"), each = q, times = q)
  phis = paste0(rep(1:q, each = q*maxLag), rep(1:q, times = q*maxLag))
  phi_pars = paste0(phi_order, phis)

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
  VARmodel = mlts_model_priors(VARmodel = VARmodel, default = T)

  # CONSTRAINTS ===============================================================

  if(q > 1 & is.null(InnoCovsZero)){
    InnoCovsZero = TRUE
  } else {
    InnoCovsZero = FALSE
  }
  if(InnoCovsZero == FALSE & q > 1 & FixInnoCovs == FALSE){
message("Note: When specifying a VAR(1) model with person-specific innovation
covariances, a latent-variable approach will be used which affords introducing
contraints on the loading parameters of the latent covariance factor(s).
(see Hamaker et al., 2018). If innovation covariances are set as a constant
or set to 0 in a subsequent step, this warning can be ignored.")
    }

  if(FixDynamics == TRUE | FixInnoVars == TRUE |
     FixInnoCovs == TRUE | !is.null(FEis0) |
     !is.null(REis0) | InnoCovsZero == T) {
    VARmodel = mlts_model_constraint(
      VARmodel = VARmodel,
      FixDynamics = FixDynamics, FixInnoVars = FixInnoVars,
      InnoCovsZero = InnoCovsZero,
      FixInnoCovs = FixInnoCovs, FEis0 = FEis0, REis0 = REis0)
  }

  # MEASUREMENT MODEL =========================================================
  if(!is.null(p)){
    VARmodel = mlts_model_measurement(
      VARmodel = VARmodel, q = q, p = p,
      btw.factor = btw.factor, btw.model = btw.model)
  }

  # BETWEEN-MODEL =============================================================
  if(!is.null(RE.pred) | !is.null(out.pred)){
    VARmodel = mlts_model_betw(VARmodel = VARmodel,
                               RE.pred =RE.pred, out.pred=out.pred,
                               out.pred.add.btw = out.pred.add.btw)
  }


  return(VARmodel)

}
