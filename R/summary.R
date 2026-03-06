#' Create a summary of a fitted model with class `mltsfit`
#'
#' @param object An object of class `mltsfit`.
#' @param priors Logical. Should priors be included in the summary?
#' Defaults to `FALSE`.
#' @param se Logical. Should the Monte Carlo Standard Error be included
#' in the summary? Defaults to `FALSE`.
#' @param prob A value between 0 and 1 to indicate the width of the credible
#' interval. Default is .95.
#' @param bpe Bayesian posterior estimate can be either "mean" (the default)
#' or the "median" of the posterior distribution.
#' @param digits Number of digits.
#' @param diag_incl_pp Logical. Should person-specific parameters be included in the
#' calculation of model convergence diagnostics (default = `FALSE`).
#' Note this is only relevant if the model was estimated with `monitor_person_pars = TRUE`.
#' @param flag_signif Add significance flags based on `prob` (default = FALSE).
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
summary.mltsfit <- function(object, priors = FALSE, se = FALSE, prob = .95,
                            bpe = c("mean"), diag_incl_pp = FALSE,
                            digits = 3, flag_signif = FALSE, ...) {

  object <- object
  model <- object$model
  class <- object$model$class
  N_obs <- object$standata$N_obs
  N_ids <- object$standata$N
  # pop_pars <- object$pop.pars.summary
  chains <- as.numeric(object$stanfit@sim$chains)
  iter <- object$stanfit@sim$iter
  thin <- object$stanfit@sim$thin
  date <- object$stanfit@date
  algorithm <- object$stanfit@stan_args[[1]]$algorithm
  bpe = ifelse(bpe == "mean", "mean", "50%")
  alpha = 1 - prob
  probs = c(alpha/2, .50, 1-alpha/2)
  prob.cols = paste0(c(100-(100-(alpha/2)*100), 100-(alpha/2)*100), "%")

  # get CIs based on prob
  # create a summary table using the monitor-function in rstan
  par_labels = mlts_param_labels(object$model)
  if(diag_incl_pp == TRUE){
    sums <- rstan::monitor(object$stanfit, probs = probs, print = FALSE)
  } else {
    sums <- rstan::monitor(
      rstan::extract(object$stanfit, pars = object$param.labels$Param_stan, permuted = FALSE, inc_warmup = TRUE),
      probs = probs, print = FALSE)
  }
  conv <- sums # backup to use for convergence checks
  sums <- round(sums[1:dim(sums)[1], 1:ncol(sums)], digits)
  sums$Param_stan = row.names(sums)
  pop_pars <- merge(par_labels, y = sums, by = "Param_stan", sort = FALSE)

  # add significance flags
  pop_pars$signif = ifelse(
    (pop_pars[,prob.cols[1]] > 0 & pop_pars[,prob.cols[2]] > 0) |
    (pop_pars[,prob.cols[1]] < 0 & pop_pars[,prob.cols[2]] < 0),
    "*", " ")

  # create prior column
  pop_pars$prior = ifelse(pop_pars$prior_type == "LKJ",
                          paste0(pop_pars$prior_type,"(",pop_pars$prior_location,")"),
                          paste0(pop_pars$prior_type,"(",pop_pars$prior_location,", ",
                                 pop_pars$prior_scale, ")"))

  # choose columns to print
  if(se == TRUE){
    colnames(pop_pars)[which(colnames(pop_pars) == "se_mean")] <- "MC.SE"
    cols = c("Param", bpe, "sd", "MC.SE", prob.cols)
  } else {
    cols = c("Param", bpe, "sd", prob.cols)
  }

  if(flag_signif == TRUE){
    cols = c(cols, "signif", "Rhat", "Bulk_ESS", "Tail_ESS")
  } else {
    cols = c(cols, "Rhat", "Bulk_ESS", "Tail_ESS")
  }

  if(priors == TRUE) {
    cols = c(cols, "prior")
  }

  infos <- mlts_model_eval(object$model)

  # print information on variables used for model estimation
  if(infos$isLatent == FALSE){
    call_inds = c(
      "Time series variables as indicated by parameter subscripts: \n",
      unlist(lapply(1:infos$q, function(x){
        paste0("  ", x, " --> ", object$standata$ts[x], "\n")
      }))
    )
  }
  if(infos$isLatent == TRUE){
    call_inds = c(
      "Time series variables as indicated by parameter subscripts: \n",
      unlist(lapply(1:infos$q, function(x){
        paste0("  ", x, " --> ",
               paste0(object$standata$ts[infos$indicators$q == x], collapse = " + "), "\n")
      }))
    )
  }

  # check for multiple groups
  if(infos$G > 1){
    call_group <- c(
      "Multiple groups detected: \n",
      unlist(lapply(1:infos$G, function(x){
        paste0("  ", x, " --> ",
               paste0(object$standata$group_lab[x]), " (N = ", object$standata$N_G[x], ") \n")
      }))
    )
  }


  # number of observations and IDs
  data_info <- paste0(
    "Data: ",
    N_obs, " observations in ", N_ids, " IDs\n"
  )

  # model convergence info
  convergence <- paste0(
    "Model convergence criteria: \n",
    "  Maximum Potential Scale Reduction Factor (PSR; Rhat): ", round(max(conv$Rhat),3), " (should be < 1.01)\n",
    "  Minimum Bulk ESS: ", min(conv$Bulk_ESS), " (should be > ", chains*100,", 100 per chain) \n",
    "  Minimum Tail ESS: ", min(conv$Tail_ESS), " (should be > ", chains*100,", 100 per chain) \n",
    "  Number of divergent transitions: ", rstan::get_num_divergent(object$stanfit),
    " (should be 0) \n"
  )



  fixef_params <- list()
  ranef_sds <- list()
  ranef_corrs <- list()
  ranef_preds <- list()
  outcomes <- list()
  outcomes_sds <- list()

  for(g in 1:infos$G){
    # get fixed effects for printing
    fixef_params[[g]] <- pop_pars[grepl("Fix", pop_pars$Type) & pop_pars$group == g, c(cols)]
    colnames(fixef_params[[g]]) <- change_colnames(fixef_params[[g]])

    # get random effects SD for printing
    ranef_sds[[g]] <- pop_pars[grepl("Random", pop_pars$Type) & pop_pars$group == g, c(cols)]
    # drop sigma_ prefix in Param
    ranef_sds[[g]][grepl("sigma_", ranef_sds$Param), "Param"] <- substr(x = ranef_sds[[g]]$Param, start = 7, stop = 30)
    colnames(ranef_sds[[g]]) <- change_colnames(ranef_sds[[g]])

    # get random effects correlation for printing
    ranef_corrs[[g]] <- pop_pars[grepl("RE correlation", pop_pars$Type) & pop_pars$group == g, c(cols)]
    # drop r_ prefix in Param
    ranef_corrs[[g]][grepl("r_", ranef_corrs[[g]]$Param), "Param"] <- gsub(
      "r_", "", ranef_corrs[[g]]$Param
    )
    colnames(ranef_corrs[[g]]) <- change_colnames(ranef_corrs[[g]])

    # get random effect predictors
    ranef_preds[[g]] <- pop_pars[grepl("RE prediction", pop_pars$Type) & pop_pars$group == g, c(cols)]
    ranef_preds[[g]][grepl(".ON.", ranef_preds[[g]]$Param), "Param"] <- gsub(
      "b_(.*).ON.(.*)", "\\1 ~ \\2", ranef_preds[[g]]$Param
    )
    colnames(ranef_preds[[g]]) <- change_colnames(ranef_preds[[g]])

    # get outcome prediction effects
    outcomes[[g]] <- pop_pars[
      grepl("Outcome prediction", pop_pars$Type) & pop_pars$group == g & !grepl("sigma_", pop_pars$Param),cols]
    outcomes[[g]]$Param <- ifelse(
      grepl(".ON.", outcomes[[g]]$Param),
      gsub("b_(.*).ON.(.*)", "\\1 ~ \\2", outcomes[[g]]$Param), ifelse(
        grepl("alpha", outcomes[[g]]$Param),
        gsub("alpha_(\\w+)", "\\1 ~ 1", outcomes[[g]]$Param),
        NA
      )
    )
    colnames(outcomes[[g]]) <- change_colnames(outcomes[[g]])

    # get outcome SDs
    outcomes_sds[[g]] <- pop_pars[grepl("Outcome prediction", pop_pars$Type) & pop_pars$group == g & grepl("sigma_", pop_pars$Param),cols]
    outcomes_sds[[g]][grepl("sigma_", outcomes_sds[[g]]$Param), "Param"] <- gsub(
    "sigma_(\\w+)", "\\Residual SD \\1", outcomes_sds[[g]]$Param
    )
    colnames(outcomes_sds[[g]]) <- change_colnames(outcomes_sds[[g]])
}

  # get measurement model parameters
  mm_pars <- pop_pars[grepl("Measurement|Item|Loading", pop_pars$Type),cols]
  colnames(mm_pars) <- change_colnames(mm_pars)



  # assemble everything
  cat(call_inds)
  if (infos$G > 1){
    cat(call_group)
  }
  cat(data_info)
  cat(convergence)

  for(g in 1:infos$G){

    if(infos$G > 1){
      group_header <- paste0("\nPosterior Summary Statistics for Group: ", object$standata$group_lab[g]," \n")
      cat(group_header)
    } else {
      cat("\nPosterior Summary Statistics")
    }

    if (nrow(fixef_params[[g]]) > 0) {
      cat("\nFixed Effects:\n")
      print(fixef_params[[g]], row.names = FALSE)
    }
    if(nrow(ranef_preds[[g]]) > 0) {
      cat("\nRandom Effects Regressed On:\n")
      print(ranef_preds[[g]], row.names = FALSE)
    }
    if(nrow(outcomes[[g]]) > 0) {
      cat("\nOutcome Prediction:\n")
      outcomes[[g]] <- rbind(outcomes[[g]], outcomes_sds[[g]])
      print(outcomes[[g]], row.names = FALSE)
    }
    if (nrow(ranef_sds[[g]]) > 0) {
      cat("\nRandom Effects SDs:\n")
      print(ranef_sds[[g]], row.names = FALSE)
    }
    if (nrow(ranef_corrs[[g]]) > 0) {
      cat("\nRandom Effects Correlations:\n")
      print(ranef_corrs[[g]], row.names = FALSE)
    }
  }

    if (nrow(mm_pars) > 0) {
      cat("\nMeasurement Model Parameters:\n")
      print(mm_pars, row.names = FALSE)
    }

    cat("\nSamples were drawn using ", algorithm, " on ", date, ".\n",
      "For each parameter, Bulk_ESS and Tail_ESS are measures of effective\n",
      "sample size, and Rhat is the potential scale reduction factor\n",
      "on split chains (at convergence, Rhat = 1).",
      sep = "")

}
