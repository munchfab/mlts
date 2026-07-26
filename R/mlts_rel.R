#' Reliability estimates for individual model parameters
#'
#' Computes reliability estimates for individual-specific parameters from a fitted
#' \code{mlts} model using approaches recommended for intensive longitudinal data.
#'
#' Reliability is conceptualized as the proportion of true individual differences
#' relative to total variance, accounting for estimation uncertainty.
#'
#' Two methods are currently implemented:
#' \itemize{
#'   \item \code{"EAP"}: Reliability based on posterior mean estimates and
#'   posterior standard deviations.
#'   \item \code{"RMU"}: Reliability based on repeated random split-half
#'   correlations of posterior draws, summarized by location and uncertainty.
#' }
#'
#' @param mltsfit An object returned by \code{mlts()}.
#' @param method Character string specifying the reliability estimation method.
#'   Must be one of \code{"EAP"} or \code{"RMU"}.
#' @param prob Numeric value between 0 and 1 specifying the probability mass
#'   covered by the equal-tailed credible interval (ETI) in the \code{"RMU"} method.
#'   Defaults to \code{0.95}.
#' @param seed Integer used in the \code{"RMU"} method.
#'
#' @return
#' If \code{method = "EAP"}:
#' \itemize{
#'   \item A named numeric vector containing reliability estimates for all
#'   random individual parameters.
#' }
#'
#' If \code{method = "RMU"}:
#' \itemize{
#'   \item A data.frame with one row per random parameter containing:
#'   \itemize{
#'     \item \code{rel_mean}: Mean reliability estimate
#'     \item \code{rel_mdn}: Median reliability estimate
#'     \item \code{rel_sd}: Standard deviation of reliability estimates
#'     \item \code{rel_ci.lb}: Lower bound of the equal-tailed credible interval
#'     \item \code{rel_ci.ub}: Upper bound of the equal-tailed credible interval
#'   }
#' }
#'
#' @details
#' The \code{"EAP"} method estimates reliability as the ratio of true score
#' variance to total variance, where total variance is decomposed into
#' between-person variability in the posterior mean estimates and the expected
#' posterior uncertainty (measurement error). Specifically, reliability is defined as
#' \deqn{
#' Rel = \frac{Var(\hat{\theta})}{Var(\hat{\theta}) + E(SE^2)}
#' }
#' where \eqn{\hat{\theta}} denotes the posterior mean (EAP) estimates of the
#' individual parameters and \eqn{SE} their posterior standard deviations. This
#' formulation follows the conceptualization of reliability as the proportion of
#' variance in observed estimates that reflects stable between-person differences,
#' rather than estimation uncertainty.
#'
#' The \code{"RMU"} method implements a resampling-based approach to quantify the
#' reliability of individual differences in model parameters. Posterior draws are
#' repeatedly split into independent subsets, and correlations between these subsets
#' are computed across individuals. These split-half correlations provide an estimate
#' of the stability (replicability) of individual differences under posterior
#' uncertainty. The resulting distribution of correlations is summarized to obtain
#' point estimates and uncertainty intervals (equal-tailed credible intervals) for
#' reliability.
#'
#' Both approaches explicitly account for uncertainty in individual parameter
#' estimation and are designed for intensive longitudinal data, where reliability
#' pertains to the stability of individual-specific model parameters rather than
#' observed variables.
#'
#' The \code{"RMU"} implementation is described in:
#' Bignardi, G., Kievit, R. A., & Bürkner, P.C. (2025),
#' \emph{A general method for estimating reliability using Bayesian Measurement Uncertainty},
#' PsyArXiv. \url{https://osf.io/preprints/psyarxiv/h54k8_v1}
#'
#' @examples
#' \dontrun{
#'
#' # Computation of individual parameter reliability for an AR(1) model:
#' ## Model specification
#' AR1 <- mlts_model(q = 1)
#'
#' ## Simulate data
#' simData <- mlts_sim(model = AR1, N = 70, TP = 70, default = TRUE)
#'
#' ## Estimate the model and store MCMC draws of individual parameters
#' fit <- mlts_fit(
#'   model = AR1, data = simData$data,
#'   id = "ID", ts = "Y1",
#'   monitor_person_pars = TRUE # !important
#' )
#'
#' # Person separation reliability
#' mlts_rel(fit, method = "EAP")
#'
#' # Relative Measurement Uncertainty (RMU)
#' mlts_rel(fit, method = "RMU")
#' }
#'
#' @export
mlts_rel <- function(mltsfit, method = c("EAP", "RMU"), prob = 0.95, seed = NULL) {
  # check that MCMC draws of individual parameter are available
  if (is.null(mltsfit$person.pars.summary$sd)) {
    stop("MCMC draws of individual parameter are not available. Refit the model with monitor_person_pars = TRUE.")
  }

  method <- match.arg(method)


  if (method == "EAP") {
    u_pars <- mltsfit$model$Param[mltsfit$model$isRandom == 1]
    split_data <- split(mltsfit$person.pars.summary, mltsfit$person.pars.summary$Param)

    reliability_results <- sapply(u_pars, function(x) {
      var_EAP <- stats::var(split_data[[x]]$mean, na.rm = TRUE)
      mean_err_var <- mean(split_data[[x]]$sd^2, na.rm = TRUE)
      rel <- var_EAP / (var_EAP + mean_err_var)
      return(rel)
    })

    return(reliability_results)
  }


  if (method == "RMU") {
    if (!is.null(seed)) {
      set.seed(seed = seed)
    }


    # new method
    draws_rstan <- rstan::extract(mltsfit$stanfit)
    b_free_array <- draws_rstan$b_free

    # 1. preparation of a list of all N*K matrices
    name_vec <- mltsfit$model$Param[mltsfit$model$Level == "Within"
                                    & mltsfit$model$isRandom == 1]
    matrix_list <- asplit(b_free_array, MARGIN = 3)
    matrix_list <- lapply(matrix_list, t) # list of all N*K matrices
    names(matrix_list) <- name_vec

    # 2. randomly sort columns and split matrices into M and W
    matrix_list <- lapply(matrix_list, function(mat) {
      mat[, sample(ncol(mat))]
    })
    split_matrix_list <- lapply(matrix_list, function(mat) {
      mid <- ncol(mat) / 2
      list(
        M = mat[, 1:mid],
        W = mat[, (mid + 1):ncol(mat)]
      )
    })

    # 3. generate vector of pearson correlations between M and W
    cor_list <- lapply(split_matrix_list, function(mat) {
      cor_vec <- sapply(1:ncol(mat[[1]]), function(i) {
        cor_val <- stats::cor(mat$M[, i], mat$W[, i], method = "pearson")
        return(cor_val)
      })
    })
    # 4. calculating reliability etc. from p
    ci_prob <- c((1 - prob) / 2, prob + (1 - prob) / 2)
    reliability_results <- do.call(rbind, lapply(cor_list, function(i) {
      ci <- stats::quantile(i, ci_prob, names = FALSE)

      data.frame(
        rel_mean = mean(i),
        rel_mdn = stats::median(i),
        rel_sd = stats::sd(i),
        rel_ci.lb = ci[1],
        rel_ci.ub = ci[2]
      )
    }))
    rownames(reliability_results) <- name_vec

    return(reliability_results)
  }
}
