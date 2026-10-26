#' Create a summary table ready for publication of a fitted model with class `mltsfit`
#'
#' @param fit An object of class `mltsfit`.
#' @param bpe Bayesian posterior estimate can be either "mean" (the default)
#' or the "median" of the posterior distribution.
#' @param ci.level A value between 0 and 1 to indicate the width of the credible
#' interval. Default is .95.
#' @param cols Character. Additional columns to be added (see `?rstan::monitor`).
#' @param digits Number of digits for fixed effects, random effects SDs, and
#' random effects correlations.
#' @param flag_signif Add significance flags based on `ci.level` (default = FALSE).
#' @param priors Add prior information (default = FALSE).
#' @param ... Additional arguments affecting the summary produced.
#'
#' @return A summary of model parameters.
#' @export
#'
#' @examples
#' \donttest{
#' # build simple vector-autoregressive mlts model for two time-series variables
#' var_model <- mlts_model(q = 2)
#'
#' # fit model with (artificial) dataset ts_data
#' fit <- mlts_fit(
#'   model = var_model,
#'   data = ts_data,
#'   ts = c("Y1", "Y2"), # time-series variables
#'   id = "ID", # identifier variable
#'   time = "time",
#'   tinterval = 1 # interval for approximation of continuous-time dynamic model,
#' )
#'
#' # inspect model summary
#' summary(fit)
#' }
mlts_tab <- function(
    fit, bpe = c("mean", "median"), ci.level = .95, cols = NULL, add_RE = TRUE,
    add_cor = TRUE, add_ranef_pred = FALSE, sort_by_RE = FALSE,
    add_cis = TRUE, add_cor_cis = TRUE, cor_diag_string = "-",
    FE_prefix = "FE__", RE_prefix = "RE__", RE_pred_prefix = "Par.ON__",
    RE_cor_prefix = "RE COR__",
    NA_string = " ", digits = 2, digits_cor_ci = 2){

  bpe <- match.arg(bpe, several.ok = F)

  alpha = 1 - ci.level
  probs = c(alpha/2, .50, 1-alpha/2)
  prob.cols = paste0(c(100-(100-(alpha/2)*100), 100-(alpha/2)*100), "%")


  if(bpe == "median"){
    bpe = "Q50"
  }

  cols = c(bpe, cols)

  if(add_cis == TRUE){
    cols = c(cols, prob.cols)
  }


  n_random = fit$standata$n_random
  n_cor_dim = ifelse(sort_by_RE == F, fit$standata$n_pars, fit$standata$n_random)

  # get results
#  sums = rstan::monitor(fit$stanfit, probs = probs, print = F)
  pars = fit$stanfit@sim$pars_oi[!(fit$stanfit@sim$pars_oi == "b_free")]
  sums = rstan::summary(fit$stanfit, pars = pars, probs = probs)$summary

  # Fixed effects ---
  params = fit$pop.pars.summary$Param[fit$pop.pars.summary$Type == "Fixed effect"]
  params.stan = fit$pop.pars.summary$Param_stan[fit$pop.pars.summary$Type == "Fixed effect"]

  rows = which(row.names(sums) %in% c(params.stan))
  FE = cbind.data.frame(params, round(sums[rows,cols],digits))
  colnames(FE)[2:(length(cols)+1)] <- paste0(FE_prefix,cols)

  # Random effects SD ---
  params = fit$pop.pars.summary$Param[fit$pop.pars.summary$Type == "Random effect SD"]
  params = gsub(params, pattern = "sigma_", replacement = "")
  params.stan = fit$pop.pars.summary$Param_stan[fit$pop.pars.summary$Type == "Random effect SD"]

  rows = which(row.names(sums) %in% c(params.stan))
  RE = cbind.data.frame(params, round(sums[rows,cols],digits))
  colnames(RE)[2:(length(cols)+1)] <- paste0(RE_prefix,cols)


  # Random effects correlation
  par_labs = params

  if(add_cor == TRUE & n_random > 1){

    # for now ...
    n_cor_dim = n_random
    cor_tab = matrix(nrow = n_cor_dim, ncol = n_cor_dim)

    for(i in 1:n_cor_dim){
      for(j in 1:n_cor_dim){
        if(i > j){
          par.name = paste0("bcorr[1,",i,",",j,"]")
          post.mean = sums[which(row.names(sums) == par.name), "mean"]
          cor_tab[i,j] = MOTE::apa(post.mean,digits,F)
        }
        if(i < j){
          if(add_cor_cis == TRUE){
            par.name = paste0("bcorr[1,",i,",",j,"]")
            CI.lb = MOTE::apa(sums[which(row.names(sums) == par.name), prob.cols[1]],digits_cor_ci, leading = F)
            CI.ub = MOTE::apa(sums[which(row.names(sums) == par.name), prob.cols[2]],digits_cor_ci, leading = F)
            cor_tab[i,j] = paste0("[",CI.lb, ", ", CI.ub,"]")
          } else {
            cor_tab[i,j] = " "
          }
        }
        if(i == j){
          cor_tab[i,j] = cor_diag_string
        }
      }
    }

    # add variable names
    colnames(cor_tab) <- paste0(RE_cor_prefix, par_labs)
    row.names(cor_tab) <- par_labs

    # add variable name for merging
    cor_tab = cbind.data.frame(params, cor_tab)
  }

  # get results of random effect predictors
  if(add_ranef_pred == TRUE & fit$standata$n_cov > 1){
    RE_pred = cbind.data.frame(params)

    # get covariates from prediction parameter labels
    params_pred = fit$pop.pars.summary$Param[fit$pop.pars.summary$Type == "RE prediction"]
    cov_labels = unique(gsub(params_pred, pattern = "b_(.*.).ON.", replacement = ""))

    for(i in 1:length(cov_labels)){
      params_pred = fit$pop.pars.summary$Param[fit$pop.pars.summary$Type == "RE prediction"]
      params_pred = params_pred[endsWith(params_pred, suffix = cov_labels[i])]
      params_pred = gsub(params_pred, pattern = ".ON(.*.)", replacement = "")
      params_pred = gsub(params_pred, pattern = "b_", replacement = "")
      params.stan = fit$pop.pars.summary$Param_stan[fit$pop.pars.summary$Type == "RE prediction"]

      rows = which(row.names(sums) %in% c(params.stan))
      RE_pred_p = cbind.data.frame(params_pred, round(sums[rows,cols],digits))
      colnames(RE_pred_p)[2:(length(cols)+1)] <- paste0(RE_pred_prefix,cov_labels[i],"_",cols)
      RE_pred = merge(x = RE_pred, y = RE_pred_p, by.x = "params", by.y = "params_pred", sort = F, all.x = T)
    }
  }


  # Combine all results in a table
  tab = cbind.data.frame(
    #  "Par" = par_labs,
    # Fixed effects
    FE)

  # Merge random effect predictors
  if(add_ranef_pred == T & fit$standata$n_cov > 1){
    tab = merge(x = tab, y = RE_pred, by = "params", sort=F, all.x = T)
  }

  # Merge Random effects
  if(add_RE == TRUE){
    tab = merge(x = tab, y = RE, by = "params", sort=F, all.x = T)
  }

  # Merge cor_tab if requested
  if(add_cor){
    tab = merge(x = tab, y = cor_tab, by = "params", sort = F, all.x = T)
  }

  # adjust names
  row.names(tab) <- NULL

  # set NAs
  if(!is.null(NA_string)){
    tab[is.na(tab)] <- NA_string
  }

  return(tab)

}

