#' Title
#'
#' @param q Integer. The number of time-varying constructs.
#' @param p Integer. For multiple-indicator models, specify a vector of length
#' `q` with the number of manifest indicators per construct.
#' @param max_lag Integer. The maximum lag of the autoregressive effect. The maximum
#' is 3.
#' @param fix_inno_vars Logical. Fix all random effect variances (except those
#' of individual traits) to zero.
#' @param fix_dynamics Logical. Fix all random effect variances (except those
#' of individual traits) to zero.
#' @param fix_inno_covs Logical. Set all innovation covariances to a constant value.
#' @param inno_covs_zero Logical. Set to TRUE to treat all innovations as independent.
#' @param fixef_zero Character. A character vector to index which fixed effects
#' should be fixed to zero (Note: this also results in removing the random effect
#' variance of the respective parameter).
#' @param ranef_zero Logical. Set to `TRUE` to treat all innovations as independent.
#' @param btw_factor Logical. If `TRUE` (the default), a common between-level factor
#' is modeled across all indicator variables. If `FALSE`, instead of a between-level
#' factor, indicator mean levels will be included as individual (random) effects drawn
#' from a joint multivariate normal distribution.
#' @param btw_model A list to indicate for which manifest indicator variables a common
#' between-level factor should be modeled (see Details for detailed instructions).
#' At this point restricted to one factor per latent construct.
#' @param ranef_pred A character vector or a named list. Include between-level covariate(s)
#' as predictor(s) of all random effects in `model` by entering a vector of unique variable
#' names. Alternatively, to include between-level covariates or differing sets of
#' between-level covariates as predictors of specific random effects, a named
#' list (using the `param`-labels in `model`) can be entered (see details).
#' Note that if a named list is provided, all names that do not match random
#' parameters in `model` will be ignored.
#' @param out_pred A character vector or a named list. Include between-level covariate(s)
#' as predictor(s) of all random effects in `model` by entering a vector of unique variable
#' names. Alternatively, to include between-level covariates or differing sets of
#' between-level covariates as predictors of specific random effects, a named
#' list (using the `param`-labels in `model`) can be entered (see details).
#' @param out_pred_add_btw tba
#' @return An object of class `data.frame`.
#' @export
#'
mlts_model <- function(q, p = NULL, max_lag = c(1,2,3),
                          btw_factor = TRUE, btw_model = NULL,
                          fix_dynamics = F, fix_inno_vars = F,
                          fix_inno_covs = T, inno_covs_zero = F,
                          fixef_zero = NULL, ranef_zero = NULL,
                          ranef_pred = NULL, out_pred=NULL, out_pred_add_btw = NULL){

  if(length(max_lag) == 3){
    max_lag = 1
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

  n_phi = (q^2)*max_lag                        # dynamic parameters
  phi_order = rep(paste0("phi(",1:max_lag,")_"), each = q, times = q)
  phis = paste0(rep(1:q, each = q*max_lag), rep(1:q, times = q*max_lag))
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
    model = rbind(FE, RE, REcors)
  }

  # ADD DEFAULT PRIORS ========================================================
  model = mlts_model_priors(model = model, default = T)

  # CONSTRAINTS ===============================================================

  if(q > 1 & is.null(inno_covs_zero)){
    inno_covs_zero = TRUE
  } else {
    inno_covs_zero = FALSE
  }
  if(inno_covs_zero == FALSE & q > 1 & fix_inno_covs == FALSE){
message("Note: When specifying a VAR(1) model with person-specific innovation
covariances, a latent-variable approach will be used which affords introducing
contraints on the loading parameters of the latent covariance factor(s).
(see Hamaker et al., 2018). If innovation covariances are set as a constant
or set to 0 in a subsequent step, this warning can be ignored.")
    }

  if(fix_dynamics == TRUE | fix_inno_vars == TRUE |
     fix_inno_covs == TRUE | !is.null(fixef_zero) |
     !is.null(ranef_zero) | inno_covs_zero == T) {
    model = mlts_model_constraint(
      model = model,
      fix_dynamics = fix_dynamics, fix_inno_vars = fix_inno_vars,
      inno_covs_zero = inno_covs_zero,
      fix_inno_covs = fix_inno_covs, fixef_zero = fixef_zero, ranef_zero = ranef_zero)
  }

  # MEASUREMENT MODEL =========================================================
  if(!is.null(p)){
    model = mlts_model_measurement(
      model = model, q = q, p = p,
      btw_factor = btw_factor, btw_model = btw_model)
  }

  # BETWEEN-MODEL =============================================================
  if(!is.null(ranef_pred) | !is.null(out_pred)){
    model = mlts_model_betw(model = model,
                               ranef_pred =ranef_pred, out_pred=out_pred,
                               out_pred_add_btw = out_pred_add_btw)
  }


  return(model)

}
