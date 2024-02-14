#' Title
#'
#' @param VARmodel data.frame. Output of VARmodel-Functions.
#'
#' @return An object of class `data.frame`.
#' @export
#'
VARmodelEval <- function(VARmodel){

  # read features of specified VARmodel

  # create additional columns
  ars = paste0("phi_", 1:9, 1:9)
  VARmodel$isAR = ifelse(VARmodel$Param %in% ars,1,0)

  # check if measurement model is entered
  isLatent = nrow(VARmodel[VARmodel$Model == "Measurement",])
  isLatent = ifelse(isLatent>0, TRUE, FALSE)

  ## Dynamic model
  fix_pars = VARmodel[VARmodel$Type == "Fix effect",]
  # extract the number of (latent constructs)
  constructs = lapply(fix_pars$Param[startsWith(fix_pars$Param_Label, "Trait")], function(x){
    substr(strsplit(x, split = c("_"))[[1]][2], start = 1, stop = 1)
  })
  q = max(as.numeric(unlist(constructs)))   # number of (latent) constructs
  n_mus = sum(startsWith(fix_pars$Param_Label,prefix = "Trait"))
  fix_pars$no = 1:nrow(fix_pars)
  fix_pars_dyn = fix_pars[fix_pars$Param_Label == "Dynamic",]
  n_dyn.pars = nrow(fix_pars_dyn)
  if(q < 10){  ### ACHTUNG AKTUELLE LÖSUNG FUNKTIONIERT NUR MIT D < 10
    fix_pars_dyn$Dout = substr(fix_pars_dyn$Param, start = 5, stop = 5)
    fix_pars_dyn$Dpred = substr(fix_pars_dyn$Param, start = 6, stop = 6)
  }
  # lagged relations between constructs
  N_pred = table(fix_pars_dyn$Dout) # number of lagged preds in each dimension
  D_pred = matrix(0, nrow = q, ncol = q, byrow = T)
  Dpos1 = c()
  Dpos2 = c()
  for(i in 1:q){
    D_pred[i,1:N_pred[i]] = as.integer(fix_pars_dyn$Dpred[fix_pars_dyn$Dout == i])
    if(i == 1){
      Dpos1[i] = n_mus+1
      Dpos2[i] = n_mus+N_pred[i]
    } else {
      Dpos1[i] = Dpos2[i-1] +1
      Dpos2[i] = Dpos2[i-1] +N_pred[i]
    }
  }


  # number of indicators per latent construct
  if (isLatent == T) {
    # separate measurement model intercepts
    alphas <- VARmodel[VARmodel$Model == "Measurement" &
                         grepl("alpha", VARmodel$Param),
                       "Param"]
    # create numeric vector with all indicators
    # and add 1 to the end (to measure number of indicators of last construct)
    ind <- c(as.numeric(gsub("(.+).(\\d)", "\\2", alphas)), 1)
    # create a difference vector
    diffs <- c(1, diff(ind))
    # extract number of indicators
    p <- ind[which(diffs <= 0) - 1]

    # extract indicator information
    ## start with within-part (which always contains all indicators)
    ind_base = extract_indicator_info(VARmodel, level = "Within", type = "Loading", incl.pos_p = T)

    ## step-wise addition: ---------------------------------------------------
    indicators = merge(
      x = ind_base, y = extract_indicator_info(VARmodel, level = "Within", type = "Measurement Error SD"), all.x = T)

    add = extract_indicator_info(VARmodel, level = "Between", type = "Item intercepts")
    if(nrow(add)>0){
      indicators = merge(x = indicators, y = add, all.x = T)
    } else {
      indicators$alpha_isFree = 0
    }
    add = extract_indicator_info(VARmodel, level = "Between", type = "Loading")
    if(nrow(add)>0){
      indicators = merge(x = indicators, y = add, all.x = T)
    } else {
      indicators$lambdaB_isFree = 0
    }
    add = extract_indicator_info(VARmodel, level = "Between", type = "Measurement Error SD")
    if(nrow(add)>0){
      indicators = merge(x = indicators, y = add, all.x = T)
      indicators$sigmaB_isFree[is.na(indicators$sigmaB_isFree)] = 0
    } else {
      indicators$sigmaB_isFree = 0
    }

    # get between-level fixed effect infos
    # base decision on the presence of indicator-specific mean values
    fixefs = VARmodel[VARmodel$Level == "Within" & startsWith(VARmodel$Param, "mu_"),]
    ind_means = extract_indicator_info(fixefs, level = "Within", type = "Fix effect")

    if(nrow(ind_means)>0){
      ind_means$mu_isFree = 1
      # add to indicator infos
      indicators = merge(indicators, y = ind_means, by = c("q","p"), all.x = T)
    }

    # add etaB label for matching
    if(!is.null(indicators$mu_isFree)){
      indicators$etaB_label = ifelse(!is.na(indicators$mu_isFree),
                                   paste0("mu_",indicators$q, ".", indicators$p),
                                   paste0("etaB_", indicators$q))
    } else {
      indicators$etaB_label = paste0("etaB_", indicators$q)
    }


    # add positions
    fixefs = VARmodel[VARmodel$Level == "Within" & startsWith(VARmodel$Param_Label, "Trait"),]
    indicators$etaB_pos = unlist(lapply(indicators$etaB_label, function(x){
      which(fixefs$Param == x)
    }))
    indicators$YB_free_pos = cumsum(indicators$sigmaB_isFree)

    # prepare infos to be passed to stan model ---------------------------------
    n_loadBfree = sum(indicators$lambdaB_isFree, na.rm = T)
    n_loadWfree = sum(indicators$lambdaW_isFree, na.rm = T)
    n_alphafree = sum(indicators$alpha_isFree, na.rm = T)
    n_sigmaBfree = sum(indicators$sigmaB_isFree, na.rm = T)
    n_sigmaWfree = sum(indicators$sigmaW_isFree, na.rm = T)

    pos_loadBfree = which(indicators$lambdaB_isFree == 1)
    pos_loadWfree = which(indicators$lambdaW_isFree == 1)
    pos_alphafree = which(indicators$alpha_isFree == 1)
    pos_sigmaBfree = which(indicators$sigmaB_isFree == 1)
    pos_sigmaWfree = which(indicators$sigmaW_isFree == 1)

    #
    n_YB_free = sum(indicators$sigmaB_isFree, na.rm = T)
    YB_free_pos = indicators$YB_free_pos
    mu_is_etaB = ifelse(indicators$sigmaB_isFree == 1 & !is.na(indicators$sigmaB_isFree), 0, 1)
    mu_etaB_pos = indicators$etaB_pos

  } else {
    p <- 1
    indicators = NA
    n_loadBfree = 0
    n_loadWfree = 0
    n_alphafree = 0
    n_sigmaWfree = 0
    n_sigmaBfree = 0
    pos_loadBfree = 0
    pos_loadWfree = 0
    pos_alphafree = 0
    pos_sigmaBfree = 0
    pos_sigmaWfree = 0
    n_YB_free = 0
    YB_free_pos = 0
    mu_is_etaB = 0
    mu_etaB_pos = 0
  }

  # which innovation variances as random effects?
  innos_rand = fix_pars[grepl(fix_pars$Param_Label, pattern="Variance"), "isRandom"]
  innos_pos = fix_pars[grepl(fix_pars$Param_Label, pattern="Variance"), "no"]
  n_innos_fix = q - sum(innos_rand)
  innos_fix_pos = cumsum(1 - innos_rand)


  n_pars = sum(VARmodel$Type == "Fix effect")
  n_random = sum(VARmodel$isRandom, na.rm = T)
  n_fixed = n_pars - n_random - n_innos_fix
  is_random = fix_pars$no[fix_pars$isRandom==1]
  is_fixed = matrix(fix_pars_dyn$no[fix_pars_dyn$isRandom==0], nrow = 1, ncol = n_fixed)
  re_pars = VARmodel[VARmodel$Type=="Fix effect" & VARmodel$isRandom==1,]
  re_pars$par_no = 1:nrow(re_pars)


  # number of innovation covariances to include
  n_inno_covs = nrow(fix_pars[grepl(fix_pars$Param_Label, pattern="Covariance"),])
  n_inno_cov_fix = n_inno_covs - sum(fix_pars[grepl(fix_pars$Param_Label, pattern="Covariance"),"isRandom"])
  inno_cov_pos = matrix(fix_pars[grepl(fix_pars$Param_Label, pattern="Covariance"), "no"],
                        nrow = 1, ncol = n_inno_covs)
  inno_cors = fix_pars[grepl(fix_pars$Param, pattern="r_zeta"),]
  n_inno_cors = nrow(inno_cors)


  # REs as OUTCOME ============================================================
  RE.PREDS = VARmodel[VARmodel$Type == "RE prediction",]

   if(nrow(RE.PREDS)>0){
    # which REs to regress on
    RE.PREDS$re_as_dv = substring(
      RE.PREDS$Param, 3,regexpr(RE.PREDS$Param,pattern = ".ON.", fixed = T)-1)
    RE.PREDS$re_preds = substring(
      RE.PREDS$Param, regexpr(RE.PREDS$Param,pattern = ".ON.", fixed = T)+4)
    RE.PREDS$re_no = sapply(RE.PREDS$re_as_dv, function(x){
      re_pars$par_no[which(re_pars$Param == x)]})
    RE.PREDS$re_pred_b_no = 1:nrow(RE.PREDS)
    re_preds_unique = unique(RE.PREDS$re_preds)
    n_cov_vars = re_preds_unique
    RE.PREDS$pred_no = sapply(RE.PREDS$re_preds, function(x){
      which(re_preds_unique == x)})

    n_cov = 1 + length(re_preds_unique)                  # add 1 for intercepts
    n_cov_bs = nrow(RE.PREDS)
    n_cov_mat = matrix(unlist(RE.PREDS[,c("pred_no", "re_no")]),
                        ncol = 2, nrow = n_cov_bs)
    n_cov_mat[,1] = n_cov_mat[,1] + 1    # shift by 1 for intercepts

   } else { # covariates
     n_cov = 1
     n_cov_bs = 0    # use a random placeholder as it will be overwritten when n_cov = 1
     n_cov_mat = matrix(1, ncol = n_cov+1, nrow = n_cov_bs)
     n_cov_vars = NA

  }


  # OUTCOME PREDICTION ========================================================
  OUT = VARmodel[VARmodel$Type=="Outcome prediction" & startsWith(VARmodel$Param, "b_"),]
  if(nrow(OUT) > 0){
    OUT$Var = substring(OUT$Param,3,regexpr(OUT$Param,pattern = ".ON.", fixed = T)-1)
    OUT$Pred = substring(OUT$Param, regexpr(OUT$Param,pattern = ".ON.", fixed = T)+4)
    OUT$Pred_Z = ifelse(!(OUT$Pred %in% fix_pars$Param), OUT$Pred, NA)
    n_z_vars = unique(OUT$Pred_Z[!is.na(OUT$Pred_Z)])
    n_z = length(n_z_vars)

    OUT$Pred_no = NA
    out_preds = unique(OUT$Pred)
    out_var = unique(OUT$Var)
    n_out = length(out_var)
    n_out_bs = matrix(ncol = 1,nrow = n_out)
    for(i in 1:n_out){
      n_out_bs[i,1] = sum(OUT$Var == out_var[i])
    }
    n_out_b_pos = matrix(0, ncol = max(n_out_bs[,1]), nrow = n_out)
    n_out_bs_sum = sum(n_out_bs[,1])
    n_out_bs_max = max(n_out_bs[,1])

    for(i in 1:n_out){
      OUT[OUT$Var == out_var[i], "out_var_no"] = i
      for(j in 1:nrow(OUT)){
        if(!is.na(OUT$Pred_Z[j])){
          OUT$Pred_no[j] = n_random + which(OUT$Pred_Z[j] == n_z_vars)
        } else{
          OUT$Pred_no[j] = re_pars[which(re_pars$Param == OUT$Pred[j]), "par_no"]
        }
      }
    }

    for(i in 1:n_out){
      n_out_b_pos[i,1:n_out_bs[i,1]] = OUT$Pred_no[OUT$Var==out_var[i]]
    }


  } else {
    n_out = 0
    n_out_bs = matrix(ncol = 1, nrow = n_out)
    n_out_bs_sum = 0
    n_out_bs_max = 0
    n_out_b_pos = matrix(nrow = n_out, ncol = 0)
    out_var = NULL
    n_z = 0
    n_z_vars = NULL
  }

  #

  # PRIORS ====================================================================
  cols = c("prior_location", "prior_scale")
  # fixed effects (intercepts)
  prior_gamma = matrix(nrow = n_random, ncol = 2)
  prior_gamma <- VARmodel[VARmodel$Type=="Fix effect" & VARmodel$isRandom == 1, cols]

  # random effect SDs
  prior_sd_R = matrix(nrow = n_random, ncol = 2)
  prior_sd_R = VARmodel[VARmodel$Type=="Random effect SD", cols]

  # random effect correlations
  prior_LKJ = unique(VARmodel$prior_location[VARmodel$Type=="RE correlation"])

  # add fix effect prior of constant innovation variance (as SD)
  prior_sigma = matrix(ncol = 2, nrow = n_innos_fix)
  if(n_innos_fix>0){
  prior_sigma = VARmodel[VARmodel$Param_Label=="Innovation Variance",cols]
  }

  # RE prediction
  prior_b_re_pred = matrix(ncol = 2, nrow = n_cov_bs)
  if(n_cov>1){
    prior_b_re_pred = RE.PREDS[,cols]
  }

  # outcome prediction
  prior_b_out = matrix(ncol = 2, nrow = n_out_bs_sum)
  prior_alpha_out = matrix(ncol = 2, nrow = n_out)
  prior_sigma_out = matrix(ncol = 2, nrow = n_out)

  if(n_out>0){
    OUT = OUT[order(OUT$out_var_no, OUT$Pred_no),]
    prior_b_out = OUT[,cols]
    for(i in 1:n_out){
      prior_alpha_out[i,1:2] = unlist(VARmodel[VARmodel$Param %in% c(paste0("alpha_",out_var[i])),cols])
      prior_sigma_out[i,1:2] = unlist(VARmodel[VARmodel$Param %in% c(paste0("sigma_",out_var[i])),cols])
    }
  }

  # END PRIOR SPECIFICATION ====================================================


  # store VARmodel with additional columns
  VARmodelext = VARmodel

  # create list to return
  VARmodelinfos = rstan::nlist(
    VARmodelext,
    q,
    p,
    n_mus,
    n_pars,
    n_random,
    n_fixed,
    is_random,
    re_pars,
    is_fixed,
    fix_pars,
    fix_pars_dyn,
    n_dyn.pars,
    N_pred,
    D_pred,
    Dpos1,
    Dpos2,
    innos_rand,
    innos_pos,
    n_innos_fix,
    innos_fix_pos,
    n_inno_covs,
    n_inno_cov_fix,
    inno_cov_pos,
    inno_cors,
    n_inno_cors,

    # REs as outcomes
    RE.PREDS, n_cov, n_cov_bs, n_cov_mat, n_cov_vars,

    # outcome model information
    OUT, n_out, out_var, n_out_bs, n_out_bs_sum, n_out_bs_max, n_out_b_pos,
    n_z, n_z_vars,

    # measurement model parameters
    isLatent, indicators, n_loadBfree, n_loadWfree, n_alphafree, n_sigmaWfree, n_sigmaBfree,
    pos_loadBfree, pos_loadWfree, pos_alphafree, pos_sigmaBfree, pos_sigmaWfree,
    n_YB_free, YB_free_pos, mu_is_etaB, mu_etaB_pos,

    # priors
    prior_gamma, prior_sd_R, prior_LKJ, prior_sigma, prior_b_re_pred,
    prior_b_out, prior_alpha_out, prior_sigma_out
    )

  return(VARmodelinfos)
}


