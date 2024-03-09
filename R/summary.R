#' Create a summary of a fitted model with class `mltsfit`
#'
#' @param object An object of class `mltsfit`.
#' @param priors Logical. Should priors be included in the summary?
#' Defaults to `FALSE`.
#' @param se Logical. Should the Monte Carlo Standard Error be included
#' in the summary? Defaults to `FALSE`.
#' @param prob A value between 0 and 1 to indicate the width of the credible
#' interval. Default is .95.
#' @param ... Additional arguments affecting the summary produced.
#'
#' @return
#' @export
#'
#' @examples
summary.mltsfit <- function(object, priors = FALSE, se = FALSE, prob = .95,
                            ...) {

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
  call_end <- ")\n"

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
        grepl("ln.sigma", infos$fix_pars$Param), "isRandom"
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

  # get fixed effects for printing
  fixef_params <- pop_pars[grepl("Fix", pop_pars$Type), c(
    "Param", "mean", "sd", "2.5%", "97.5%", "Rhat", "Bulk_ESS", "Tail_ESS"
  )]
  colnames(fixef_params)[1:3] <- c("", "Estimate", "Post.SD")

  # get random effects for printing
  ranef_params <- pop_pars[grepl("Random", pop_pars$Type), c(
    "Param", "mean", "sd", "2.5%", "97.5%", "Rhat", "Bulk_ESS", "Tail_ESS"
  )]
  # drop sigma_ prefix in Param
  ranef_params[grepl("sigma_", ranef_params$Param), "Param"] <- gsub(
    "sigma_", "", ranef_params$Param
  )
  colnames(ranef_params)[1:3] <- c("", "Estimate", "Post.SD")

  # assemble everything
  cat(model_call)
  if (nrow(fixef_params) > 0) {
    cat("\nFixed Effects:\n")
    print(fixef_params, row.names = FALSE)
  }
  if (nrow(ranef_params) > 0) {
    cat("\nRandom Effects SD:\n")
    print(ranef_params, row.names = FALSE)
  }


}
