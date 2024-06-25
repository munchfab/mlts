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
                            bpe = c("mean", "median"),
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
  sums <- rstan::monitor(object$stanfit, probs = probs, print = FALSE)
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

  # assemble model call to paste in summary
  call_start <- paste0(
    "Call:\n",
    "mlts_model(q = ", infos$q,
    sep = ""
  )
  call_end <- ")\n"

  # for latent variabels: print number of indicators and latent variables
  call_p <- ifelse(
    infos$isLatent == TRUE,
    paste0(
      ", p = ", ifelse(
        length(infos$p) > 1,
        paste0("c(", paste(infos$p, collapse = ", "), ")"),
        infos$p
      )
    ), ""
  )

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

  # if (infos$isLatent == TRUE) {
  #   latents <- pop_pars[
  #     startsWith(pop_pars$Param, "eta") | startsWith(pop_pars$Param, "mu"),
  #     "Param"
  #   ]
  #   n_latents <- length(latents)
  #   indicators <- rownames(object$standata$y)
  #   # n_indicators <- length(indicators)
  #   lat_ind <- data.frame(
  #     latents = rep(latents, times = infos$p),
  #     by = "=~",
  #     indicators = indicators
  #   )
  #   mm <- aggregate(
  #     indicators ~ latents + by, data = lat_ind,
  #     FUN = function(x) {paste(x, collapse = " + ")}
  #   )
  #   mm_string <- do.call(paste, mm)
  #   call_latents <- c(
  #     paste0("Latent Variables: ", n_latents, "\n"),
  #     paste0("  ", mm_string, collapse = "\n "),
  #     "\n"
  #   )
  # } else {
  #   latents <- pop_pars[
  #     startsWith(pop_pars$Param, "mu"),
  #     "Param"
  #   ]
  #   n_latents <- length(latents)
  #   indicators <- rownames(object$standata$y)
  #   # n_indicators <- length(indicators)
  #   lat_ind <- data.frame(
  #     latents = rep(latents, times = infos$p),
  #     by = "=~",
  #     indicators = indicators
  #   )
  #   mm <- aggregate(
  #     indicators ~ latents + by, data = lat_ind,
  #     FUN = function(x) {paste(x, collapse = " + ")}
  #   )
  #   mm_string <- do.call(paste, mm)
  #   call_latents <- c(
  #     paste0("Latent Variables: ", n_latents, "\n"),
  #     paste0("  ", mm_string, collapse = "\n "),
  #     "\n"
  #   )
  # }

  call_maxlag <- paste0(", max_lag = ", infos$maxLag)
  call_fix_dynamics <- ifelse(
    all(infos$fix_pars_dyn$isRandom == 0),
    paste0(", fix_dynamics = TRUE"),
    ""
  )
  call_fix_inno_vars <- ifelse(
    all(
      infos$fix_pars[
        grepl("ln.sigma", infos$fix_pars$Param), "isRandom"
      ] == 0
    ),
    paste0(", fix_inno_vars = TRUE"),
    ""
  )
  call_ranef_zero <- ifelse(
    # determine if there is exactly one random effect == 0
    # apart from residual correlations (which are fixed by default)
    sum(
      infos$fix_pars[!grepl("r.zeta", infos$fix_pars$Param), "isRandom"] == 0
    ) == 1,
    paste0(", ranef_zero = \"",
           infos$fix_pars[!grepl("r.zeta", infos$fix_pars$Param) &
                            infos$fix_pars$isRandom == 0, "Param"],
           "\""),
    ifelse(
      # determine if there are more than one random effects == 0
      # apart from residual correlations (which are fixed by default)
      sum(
        infos$fix_pars[!grepl("r.zeta", infos$fix_pars$Param), "isRandom"] == 0
      ) > 1,
      paste0(
        ", ranef_zero = c(",
        paste0(
          "\"", infos$fix_pars[!grepl("r.zeta", infos$fix_pars$Param) &
                                 infos$fix_pars$isRandom == 0, "Param"], "\"",
          collapse = ", "
        ),
        ")"
      ), ""
    )
  )
  call_ranef_pred <- ifelse(
    length(infos$n_cov_vars) > 1,
    paste0(
      ", ranef_pred = ", ifelse(
        length(infos$n_cov_vars) == 1,
        paste0("\"", infos$n_cov_vars, "\""),
        paste0("c(", paste0(paste0("\"", infos$n_cov_vars, "\""), collapse = ", "),
               ")")
      )
    ), ""
  )
  call_out_pred <- ifelse(
    length(infos$out_var) > 0,
    paste0(
      ", out_pred = ", ifelse(
        length(infos$out_var) == 1,
        paste0("\"", infos$out_var, "\""),
        paste0("c(", paste0(paste0("\"", infos$out_var, "\""), collapse = ", "),
               ")")
      )
    ), ""
  )
  # determine fixed effects fixed to zero
  suppressMessages({
  sat_model <- mlts_model(
    q = infos$q,
    p = if (all(infos$p == 1)) {NULL} else {infos$p},
    max_lag = infos$maxLag,
    ranef_pred = if (length(infos$n_cov_vars) > 1) {infos$n_cov_vars} else {NULL},
    out_pred = if (infos$n_out > 0) {infos$n_out} else {NULL},
    fix_inno_covs =  if (infos$n_inno_cov_fix > 0) {TRUE} else {FALSE},
    inno_covs_zero =  if(infos$n_inno_covs == 0) {TRUE} else {FALSE},
    inno_covs_dir = if(is.na(infos$inno_cov_dir)) {NULL} else {infos$inno_cov_dir},
  )}
  )

  sat_model_fixed <- sat_model[grepl("Fix", sat_model$Type), "Param"]
  model_fixed <- model[grepl("Fix", model$Type), "Param"]
  fe_zero <- setdiff(
    union(sat_model_fixed, model_fixed), intersect(sat_model_fixed, model_fixed)
  )
  call_fixef_zero <- ifelse(
    length(fe_zero) > 0,
    paste0(
      ", fixef_zero = ", ifelse(
        length(fe_zero) == 1,
        paste0("\"", fe_zero, "\""),
        paste0("c(", paste0(paste0("\"", fe_zero, "\""), collapse = ", "),
               ")")
      )
    ), ""
  )


  # concatenate and paste to summary output
  model_call <- paste(
    call_start, call_p, call_maxlag,
    # maybe not necessary?
    call_fix_dynamics, call_fix_inno_vars,
    call_fixef_zero, call_ranef_zero,
    call_ranef_pred, call_out_pred,
    call_end,
    sep = ""
  )

  # number of observations and IDs
  data_info <- paste0(
    "Data: ",
    N_obs, " observations in ", N_ids, " IDs\n"
  )

  # model convergence info
  conv = rstan::monitor(object$stanfit, print = FALSE, probs = prob)
  convergence <- paste0(
    "Model convergence criteria: \n",
    "  Maximum Potential Scale Reduction Factor (PSR; Rhat): ", round(max(conv$Rhat),3), " (should be < 1.01)\n",
    "  Minimum Bulk ESS: ", min(conv$Bulk_ESS), " (should be > ", chains*100,", 100 per chain) \n",
    "  Minimum Tail ESS: ", min(conv$Tail_ESS), " (should be > ", chains*100,", 100 per chain) \n",
    "  Number of divergent transitions: ", rstan::get_num_divergent(object$stanfit),
    " (should be 0) \n"
  )

  # Type is missing pop_pars if model is latent
  # disable this workaround when fixed
  # if (infos$isLatent == TRUE) {
  #   pop_pars <- merge(pop_pars, object$model[, c("Type", "Param")], by = "Param")
  # }


  # get fixed effects for printing
  fixef_params <- pop_pars[grepl("Fix", pop_pars$Type), c(cols)]
  # colnames(fixef_params)[1:3] <- c("", "Estimate", "Post.SD")
  colnames(fixef_params) <- change_colnames(fixef_params)

  # get random effects SD for printing
  ranef_sds <- pop_pars[grepl("Random", pop_pars$Type), c(cols)]
  # drop sigma_ prefix in Param
  ranef_sds[grepl("sigma_", ranef_sds$Param), "Param"] <- substr(x = ranef_sds$Param, start = 7, stop = 30)
  # colnames(ranef_sds)[1:3] <- c("", "Estimate", "Post.SD")
  colnames(ranef_sds) <- change_colnames(ranef_sds)

  # get random effects correlation for printing
  ranef_corrs <- pop_pars[grepl("RE correlation", pop_pars$Type), c(cols)]
  # drop r_ prefix in Param
  ranef_corrs[grepl("r_", ranef_corrs$Param), "Param"] <- gsub(
    "r_", "", ranef_corrs$Param
  )
  # colnames(ranef_corrs)[1:3] <- c("", "Estimate", "Post.SD")
  colnames(ranef_corrs) <- change_colnames(ranef_corrs)

  # get random effect predictors
  ranef_preds <- pop_pars[grepl("RE prediction", pop_pars$Type), c(cols)]
  ranef_preds[grepl(".ON.", ranef_preds$Param), "Param"] <- gsub(
    "b_(.*).ON.(.*)", "\\1 ~ \\2", ranef_preds$Param
  )
  # colnames(ranef_preds)[1:3] <- c("", "Estimate", "Post.SD")
  colnames(ranef_preds) <- change_colnames(ranef_preds)

  # get outcome prediction effects
  outcomes <- pop_pars[
    grepl("Outcome prediction", pop_pars$Type) & !grepl("sigma_", pop_pars$Param),cols]
  outcomes$Param <- ifelse(
    grepl(".ON.", outcomes$Param),
    gsub("b_(.*).ON.(.*)", "\\1 ~ \\2", outcomes$Param), ifelse(
      grepl("alpha", outcomes$Param),
      gsub("alpha_(\\w+)", "\\1 ~ 1", outcomes$Param),
      NA
    )
  )
  # colnames(outcomes)[1:3] <- c("", "Estimate", "Post.SD")
  colnames(outcomes) <- change_colnames(outcomes)

  # get outcome SDs
  outcomes_sds <- pop_pars[grepl("Outcome prediction", pop_pars$Type) & grepl("sigma_", pop_pars$Param),cols]
  outcomes_sds[grepl("sigma_", outcomes_sds$Param), "Param"] <- gsub(
    "sigma_(\\w+)", "\\Residual SD \\1", outcomes_sds$Param
  )
  # colnames(outcomes_sds)[1:3] <- c("", "Estimate", "Post.SD")
  colnames(outcomes_sds) <- change_colnames(outcomes_sds)

  # get measurement model parameters
  mm_pars <- pop_pars[grepl("Measurement|Item|Loading", pop_pars$Type),cols]
  # colnames(mm_pars)[1:3] <- c("", "Estimate", "Post.SD")
  colnames(mm_pars) <- change_colnames(mm_pars)



  # assemble everything
  cat(model_call)
  # cat(call_latents)
  cat(call_inds)
  cat(data_info)
  cat(convergence)
  if (nrow(fixef_params) > 0) {
    cat("\nFixed Effects:\n")
    print(fixef_params, row.names = FALSE)
  }
  if (nrow(ranef_sds) > 0) {
    cat("\nRandom Effects SDs:\n")
    print(ranef_sds, row.names = FALSE)
  }
  if (nrow(ranef_corrs) > 0) {
    cat("\nRandom Effects Correlations:\n")
    print(ranef_corrs, row.names = FALSE)
  }
  if(nrow(outcomes) > 0) {
    cat("\nOutcome Prediction:\n")
    outcomes <- rbind(outcomes, outcomes_sds)
    print(outcomes, row.names = FALSE)
  }
  if(nrow(ranef_preds) > 0) {
    cat("\nRandom Effects Regressed On:\n")
    print(ranef_preds, row.names = FALSE)
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
