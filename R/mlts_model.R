#' Build a multilevel latent time series model
#'
#' @param class Character. Indicating the model type to be specified. For now
#' restricted to `VAR`, the default. Future package releases might include additional
#' model types.
#' @param q Integer. The number of time-varying constructs.
#' @param p Integer. For multiple-indicator models, specify a vector of length
#' `q` with the number of manifest indicators per construct. If all constructs are
#' measured with the same number of indicators, a single value is sufficient.
#' @param max_lag Integer. The maximum lag of the autoregressive effect to be
#' included in the model. The maximum is 3. Defaults to 1.
#' @param fix_dynamics Logical. Fix all random effect variances of autoregressive and
#' cross-lagged effects to zero (constraining parameters to be equal across clusters).
#' @param fix_inno_vars Logical. Fix all random effect variances of innovation variances
#' to zero (constraining parameters to be equal across clusters).
#' @param fix_inno_covs Logical. Fix all random effect variances of innovation covariances
#' to zero (constraining parameters to be equal across clusters).
#' @param inno_covs_zero Logical. Set to `TRUE` to treat all innovations as independent.
#' @param inno_covs_dir For bivariate VAR models with person-specific innovation covariances,
#' a latent variable approach is applied (for a detailed description, see Hamaker et al., 2018).
#' by specifying an additional factor that loads onto the contemporaneous innovations of both constructs,
#' capturing the shared variance of innovations, that is not predicted by the previous time points.
#' The loading parameters of this latent factor, however, have to be restricted in accordance with
#' researchers assumptions about the sign of the association between innovations across construct.
#' Hence, if innovations at time $t$ are assumed to be positively correlated across clusters, set the
#' argument to `pos`, or `neg` respectively.
#' @param fixef_zero Character. A character vector to index which fixed effects
#' (referring to the parameter labels in `model$Param`) should be constrained to zero
#' (Note: this also results in removing the random effect variance of the respective parameter).
#' @param ranef_zero Character. A character vector to index which random effect variances
#' (referring to the parameter labels in `model$Param`) should be constrained to zero.
#' @param btw_factor Logical. If `TRUE` (the default), a common between-level factor
#' is modeled across all indicator variables per construct `q`. If `FALSE`, instead of a between-level
#' factor, indicator mean levels will be included as individual (random) effects drawn
#' from a joint multivariate normal distribution.
#' @param btw_model A list to indicate for which manifest indicator variables a common
#' between-level factor should be modeled (see Details for detailed instructions).
#' At this point restricted to one factor per latent construct.
#' @param ranef_pred A character vector or a named list. Include between-level covariate(s)
#' as predictor(s) of all random effects in `model` by entering a vector of unique variable
#' names. Alternatively, to include between-level covariates or differing sets of
#' between-level covariates as predictors of specific random effects, a named
#' list (using the labels in `model$Param`) can be entered (see examples).
#' Note that if a named list is provided, all names that do not match random
#' parameters in `model` will be ignored. Note that variables entered in `ranef_pred` will
#' be grand-mean centered by default when fitting the model with `mlts_fit`.
#' @param out_pred A character vector or a named list. Include between-level outcome(s)
#' to be regressed on all random effects in `model` by entering a vector of unique variable
#' names. Alternatively, to include multiple between-level outcomes regressed differing sets of
#' specific random effects, a named list (using the labels in `model$Param`) can be entered
#' (see examples). Note that if a named list is provided, all character strings in the vector of each list
#' (with independent variables) element that do not match random effect parameter names
#' in `model$Param` will be treated as additional between-level predictors.
#' @param out_pred_add_btw A character vector. If `out_pred` is a character (vector), all
#' inputs will be treated as between-level covariates to be used as additional predictors of
#' all outcomes specified in `out_pred`.
#' @return An object of class `data.frame` with the following columns:
#' \item{Model}{Indicates if the parameter in the respective row is part of the structural, or
#' the measurement model (if multiple indicators per construct are provided)}
#' \item{Level}{Parameter on the between- or within-level.}
#' \item{Type}{Describes the parameter type.}
#' \item{Param}{Parameter names to be referred to in arguments of `mlts_model`.}
#' \item{Param_Label}{Parameter labels (additional option to address specific parameters).}
#' \item{isRandom}{Indicates which within-level parameters are modeled as random (1) or a constant
#' across clusters (0).}
#' \item{Constraint}{Optional. Included if multiple-indicators per construct (p > 1) are provided.
#' Constraints on measurement model parameters can be changed by overwriting the respective value
#' in `model`. Possible inputs are "free", "= 0" (for SDs of measurement error variances),
#' and "= 1" (for loading parameters).}
#' \item{prior_type}{Contains the parameters' prior distribution used in `mlts_fit` (prior classes
#' can not be changed at this point).}
#' \item{prior_location}{Location values of the parameters' prior distribution used
#' in `mlts_fit` (can be changed to any real value by overwriting the respective value in `model`).}
#' \item{prior_scale}{Scale values of the parameters' prior distribution used
#' in `mlts_fit` (can be changed to any real value by overwriting the respective value in `model`).}
#' @references
#' Hamaker, E. L., Asparouhov, T., Brose, A., Schmiedek, F., & Muthén, B. (2018).
#' At the frontiers of modeling intensive longitudinal data: Dynamic structural equation models
#' for the affective measurements from the COGITO study. *Multivariate behavioral research*, *53*(6), 820-841.
#' \doi{10.1080/00273171.2018.1446819}
#' @export
#'
#' @examples
#' \donttest{
#'  # To illustrate the general model building procedure, starting with a simple
#'  # two-level AR(1) model with person-specific individual means, AR effects,
#'  # and innovation variances (the default option when using mlts_model() and q = 1).
#'  model <- mlts_model(q = 1)
#'
#'  # All model parameters (with their labels stored in model$Param) can be inspected by calling:
#'  model
#'
#'  # Possible model extensions/restrictions:
#'  # 1. Introducing additional parameter constraints, such as fixing specific
#'  #    parameters to a constant value by setting the respective random effect
#'  #    variances to zero, such as e.g. (log) innovation variances
#'  model <- mlts_model(q = 1, ranef_zero = "ln.sigma2_1")
#'  #    Note that setting the argument `fix_inno_vars` to `TRUE` provides
#'  #    a shortcut to fixing the innovation variances of all constructs
#'  #    (if q >= 1) to a constant.
#'
#'  # 2. Including a multiple indicator model, where the construct is measured by
#'  #    multiple indicators (here, p = 3 indicators)
#'  model <- mlts_model(
#'           q = 1, # the number of time-varying constructs
#'           p = 3, # the number of manifest indicators
#'           # assuming a common between-level factor (the default)
#'           btw_factor = TRUE
#'         )
#'
#'  # 3. Incorporating between-level variables. For example, inclusion of
#'  #    an additional between-level variable ("cov1") as predictor of all
#'  #    (ranef_pred = "cov1") or a specific set of random effects
#'  #    (ranef_pred = list("phi(1)_11") = "cov1"), an external outcome (e.g., "out1")
#'  #    to be predicted by all (out_pred = "out1") or specific random effects
#'  #    (out_pred = list("out1" = c("etaB_1", "phi(1)_11")), using the latent
#'  #    between-level factor trait scores (etaB_1) and individual first-order
#'  #    autoregressive effects (phi(1)_11) as joint predictors of outcome "out1".
#'  model <- mlts_model(
#'             q = 1,
#'             p = 3,
#'             fix_inno_vars = TRUE,
#'             ranef_pred = "cov1",
#'             out_pred = list("out1" = c("etaB_1", "phi(1)_11"))
#'            )
#'  #    Note that the names of the random effect parameters must match the
#'  #    parameter labels provided in model$Param, the result of the
#'  #    mlts_model()-functions.
#'
#' }

mlts_model <- function(class = c("VAR"), q, p = NULL, max_lag = c(1,2,3),
                          btw_factor = TRUE, btw_model = NULL,
                          fix_dynamics = FALSE, fix_inno_vars = FALSE,
                          fix_inno_covs = TRUE, inno_covs_zero = FALSE,
                          inno_covs_dir = NULL,
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


  if(q == 2 & inno_covs_zero == FALSE & fix_inno_covs == FALSE){
    if(is.null(inno_covs_dir)){
      stop(
        "For a bivariate VAR model with person-specific innovation covariances, ",
        "a latent variable appraoch will be used. This affords putting a restriction ",
        "on the loading parameters of the latent innovation covariance factor,",
        "specifying the association of innovations as either positive (inno_covs_dir == 'pos') ",
        "or negative (inno_covs_dir = 'neg').")
    }
    if(inno_covs_dir == "pos"){

    } else if(inno_covs_dir == "neg"){

    } else {
      stop(
        "Inputs of `inno_covs_dir` should be one of 'pos' or 'neg'"
      )
    }
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
      "Type" = "Fixed effect",
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
      "Type" = "Fixed effect",
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

  # ---
  if(q == 2 & inno_covs_zero == FALSE & fix_inno_covs == FALSE){
    model$Constraint[model$Type == "Fixed effect" & grepl(pattern = "Covariance", model$Param_Label)] = inno_covs_dir
  }


  # ADD DEFAULT PRIORS ========================================================
  model = mlts_model_priors(model = model, default = TRUE)

  # CONSTRAINTS ===============================================================

  if(q > 1 & inno_covs_zero == TRUE){
    inno_covs_zero = TRUE
  } else {
    inno_covs_zero = FALSE
  }
  if(inno_covs_zero == FALSE & q > 2 & fix_inno_covs == FALSE){
stop("Note: At this point, specifying a VAR(1) model with person-specific innovation covariances",
     "is only possible for `q == 2`.")
    }

  if(fix_dynamics == TRUE | fix_inno_vars == TRUE |
     fix_inno_covs == TRUE | !is.null(fixef_zero) |
     !is.null(ranef_zero) | inno_covs_zero == TRUE) {
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

  # we should think about something like this:
  # adding the class to model, or ideally as a specific class of data frame
  # referring to the specific models (to come), however, just adding a class to the
  # data frame resulted in turning it into a list object. The latter code lead to problems
  # in the summary.
  #return(cbind("class" = class, model))
}

