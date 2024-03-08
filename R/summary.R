summary.mltsfit <- function(object, priors = FALSE, se = FALSE, prob = .95) {

  object <- object
  pop_pars <- object$pop.pars.summary
  chains <- object$stanfit@sim$chains
  iter <- object$stanfit@sim$iter
  thin <- object$stanfit@sim$thin
  date <- object$stanfit@date
  algorithm <- object$stanfit@stan_args[[1]]$algorithm

  infos <- mlts_model_eval(object$model)

  # assemble model call to paste in summary
  call_start <- paste0(
    "Call:\n",
    "mlts_model(q = ", infos$q,
    sep = ""
  )
  call_end <- ")"

  call_p <- ifelse(
    infos$isLatent == TRUE,
    paste0(
      ", p = ", ifelse(
        length(infos$p) == 1,
        infos$p,
        cat("c(", cat(infos$p, sep = ", "), ")", sep = "")
      )
    ), ""
  )
  call_maxlag <- paste0(", max_lag = ", infos$maxLag)
  call_fix_dynamics <- ifelse(
    all(infos$fix_pars_dyn$isRandom == 0),
    paste0(", fix_dynamics = TRUE"),
    ""
  )
  call_fix_inno_vars <- ifelse(
    all(
      infos$fix_pars[
        infos$fix_pars$Param_Label == "Innovation Variance", "isRandom"
      ] == 0
    ),
    paste0(", fix_inno_vars = TRUE"),
    ""
  )

  # concatenate and paste to summary output
  model_call <- paste(
    call_start, call_p, call_maxlag,
    # maybe not necessary?
    call_fix_dynamics, call_fix_inno_vars,
    call_end,
    sep = ""
  )

  cat(model_call)

}
