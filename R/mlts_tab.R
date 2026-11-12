
#' Create a publication-ready summary table of a fitted `mltsfit` model
#'
#' Create a summary table containing posterior estimates and credible
#' intervals for fixed effects, random effect standard deviations,
#' random effect correlations, and optionally predictors of random effects.
#'
#' @param fit An object of class `mltsfit`.
#' @param bpe Character. The Bayesian posterior estimate to report.
#'   Either `"mean"` (the default) or `"median"`.
#' @param ci.level Numeric. The desired width of the credible interval,
#'   specified as a value between 0 and 1. Defaults to `.95`.
#' @param cols Character vector. Additional summary statistics to include
#'   in the table. These correspond to the column names returned by
#'   `rstan::summary()`, such as `"sd"`, `"se_mean"`, `"n_eff"`, and
#'   `"Rhat"`. The selected posterior estimate and credible interval
#'   columns are added automatically.
#' @param add_RE Logical. Whether to include random effect standard
#'   deviations in the summary table. Defaults to `TRUE`.
#' @param add_cor Logical. Whether to include the correlation matrix
#'   of the random effects. Defaults to `TRUE`.
#' @param add_ranef_pred Logical. Whether to include predictors of
#'   random effects in the summary table. Defaults to `FALSE`.
#' @param add_cis Logical. Whether to include posterior credible
#'   intervals for the reported parameters. Defaults to `TRUE`.
#' @param add_cor_cis Logical. Whether to report credible intervals
#'   for random effect correlations. Defaults to `TRUE`.
#'   If `FALSE`, the upper triangular part of the correlation matrix
#'   is left blank.
#' @param cor_diag_string Character. String used to represent the
#'   diagonal elements of the random effect correlation matrix.
#'   Defaults to `"-"`.
#' @param FE_prefix Character. Prefix added to the column names
#'   containing fixed effect estimates and summary statistics.
#'   Defaults to `"FE__"`.
#' @param RE_prefix Character. Prefix added to the column names
#'   containing random effect standard deviations and summary
#'   statistics. Defaults to `"RE__"`.
#' @param RE_pred_prefix Character. Prefix added to the column names
#'   containing estimates of predictors of random effects.
#'   Defaults to `"Par.ON__"`.
#' @param RE_cor_prefix Character. Prefix added to the column names
#'   of the random effect correlation matrix.
#'   Defaults to `"RE COR__"`.
#' @param NA_string Character or `NULL`. String used to replace
#'   missing values in the resulting table. If `NULL`, missing values
#'   are retained as `NA`. Defaults to `" "`.
#' @param digits Integer. Number of decimal places used for fixed
#'   effects, random effect standard deviations, and random effect
#'   correlations. Defaults to `2`.
#' @param digits_cor_ci Integer. Number of decimal places used for
#'   the lower and upper bounds of random effect correlation
#'   credible intervals. Defaults to `2`.
#' @param out_as_list Logical. Whether to return the summary tables
#'   as a list. If `TRUE`, a list containing one table per group
#'   is returned. If `FALSE`, the group-specific tables are
#'   combined into a single data frame. Defaults to `FALSE`.
#'
#' @return A data frame or a list of data frames containing posterior
#'   summaries of the model parameters. Depending on the selected
#'   arguments, the output can include fixed effects, random effect
#'   standard deviations, random effect correlations, and predictors
#'   of random effects. For models with multiple groups, the group
#'   identifier is included in the output.
#'
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
#' )
#'
#' # prepare the results table
#' mlts_tab(fit)
#' }
mlts_tab <- function(
    fit, bpe = c("mean", "median"), ci.level = .95, cols = NULL, add_RE = TRUE,
    add_cor = TRUE, add_ranef_pred = FALSE,
    add_cis = TRUE, add_cor_cis = TRUE, cor_diag_string = "-",
    FE_prefix = "FE__", RE_prefix = "RE__", RE_pred_prefix = "Par.ON__",
    RE_cor_prefix = "RE COR__",
    NA_string = " ", digits = 2, digits_cor_ci = 2, out_as_list = TRUE){

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


  tab_list = list()

  n_random = fit$standata$n_random
  n_group  = fit$standata$G
  if(n_group == 1){
    fit$pop.pars.summary$group = 1
  }

  # get results
  pars = fit$stanfit@sim$pars_oi[!(fit$stanfit@sim$pars_oi == "b_free")]
  sums = rstan::summary(fit$stanfit, pars = pars, probs = probs)$summary


  for(gg in 1:n_group){

    # Fixed effects ---
    params      = fit$pop.pars.summary$Param[     fit$pop.pars.summary$Type == "Fixed effect" & fit$pop.pars.summary$group == gg]
    params.stan = fit$pop.pars.summary$Param_stan[fit$pop.pars.summary$Type == "Fixed effect" & fit$pop.pars.summary$group == gg]

    rows = which(row.names(sums) %in% c(params.stan))
    FE = cbind.data.frame(params, round(sums[rows,cols],digits))
    colnames(FE)[2:(length(cols)+1)] <- paste0(FE_prefix,cols)

    # Random effects SD ---
    params = fit$pop.pars.summary$Param[fit$pop.pars.summary$Type == "Random effect SD" & fit$pop.pars.summary$group == gg]
    params = gsub(params, pattern = "sigma_", replacement = "")
    params.stan = fit$pop.pars.summary$Param_stan[fit$pop.pars.summary$Type == "Random effect SD" & fit$pop.pars.summary$group == gg]

    rows = which(row.names(sums) %in% c(params.stan))
    RE = cbind.data.frame(params, round(sums[rows,cols],digits))
    colnames(RE)[2:(length(cols)+1)] <- paste0(RE_prefix,cols)


    # Random effects correlation
    par_labs = params

    if(add_cor == TRUE & n_random > 1){

      n_cor_dim = n_random
      cor_tab = matrix(nrow = n_cor_dim, ncol = n_cor_dim)

      for(i in 1:n_cor_dim){
        for(j in 1:n_cor_dim){
          if(i > j){
            par.name = paste0("bcorr[",gg,",",i,",",j,"]")
            post.mean = sums[which(row.names(sums) == par.name), "mean"]
            cor_tab[i,j] = MOTE::apa(post.mean,digits,F)
          }
          if(i < j){
            if(add_cor_cis == TRUE){
              par.name = paste0("bcorr[",gg,",",i,",",j,"]")
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
      params_pred = fit$pop.pars.summary$Param[fit$pop.pars.summary$Type == "RE prediction" & fit$pop.pars.summary$group == gg]
      cov_labels = unique(gsub(params_pred, pattern = "b_(.*.).ON.", replacement = ""))

      for(i in 1:length(cov_labels)){
        params_pred = fit$pop.pars.summary$Param[fit$pop.pars.summary$Type == "RE prediction" & fit$pop.pars.summary$group == gg]
        params_pred = params_pred[endsWith(params_pred, suffix = cov_labels[i])]
        params_pred = gsub(params_pred, pattern = ".ON(.*.)", replacement = "")
        params_pred = gsub(params_pred, pattern = "b_", replacement = "")
        params.stan = fit$pop.pars.summary$Param_stan[fit$pop.pars.summary$Type == "RE prediction" & fit$pop.pars.summary$group == gg]

        rows = which(row.names(sums) %in% c(params.stan))
        RE_pred_p = cbind.data.frame(params_pred, round(sums[rows,cols],digits))
        colnames(RE_pred_p)[2:(length(cols)+1)] <- paste0(RE_pred_prefix,cov_labels[i],"_",cols)
        RE_pred = merge(x = RE_pred, y = RE_pred_p, by.x = "params", by.y = "params_pred", sort = TRUE, all.x = T)
      }
    }


  # Combine all results in a table
  tab = cbind.data.frame(
    #  "Par" = par_labs,
    # Fixed effects
    FE)

  # Merge random effect predictors
  if(add_ranef_pred == TRUE & fit$standata$n_cov > 1){
    tab = merge(x = tab, y = RE_pred, by = "params", all.x = TRUE)
  }

  # Merge Random effects
  if(add_RE == TRUE){
    tab = merge(x = tab, y = RE, by = "params", all.x = TRUE)
  }

  # Merge cor_tab if requested
  if(add_cor == TRUE & n_random > 1){
    tab = merge(x = tab, y = cor_tab, by = "params", all.x = TRUE)
  }

  # adjust names
  row.names(tab) <- NULL

  # set NAs
  if(!is.null(NA_string)){
    tab[is.na(tab)] <- NA_string
  }


  tab_list[[gg]] = tab

  } # end loop over groups

  if( n_group > 1 ){
    group_labels = fit$standata$group_lab

    tab_list = lapply(1:length(tab_list), function(x){
        cbind("Group" = group_labels[x], tab_list[[x]])
      })
    }



  if(out_as_list == FALSE){
    tab_out = do.call(rbind, tab_list)
  } else {
    tab_out = tab_list
  }



  return(tab_out)

}

