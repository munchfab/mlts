#' Posterior Predictive Checks for Multilevel Latent Time Series Models
#'
#' @description
#' {Developemental} This function plots posterior predictive distributions of one or multiple fitted \code{mlts.fit} models.
#' Simulated data from posterior draws are compared to the observed data to visually assess model fit.
#'
#' @param fit A fitted model object of class \code{mlts.fit}. Only used if \code{fit_list} is \code{NULL}.
#' @param fit_list An optional list of fitted \code{mlts.fit} objects for model comparison. If provided, \code{fit} is ignored.
#' @param ts Optional vector of variable names to include in the plot.
#' @param y_reps Optional. A list of posterior predictive samples (as returned by \code{\link{mlts_posterior_sample}}) for a single model.
#' If \code{NULL}, samples are generated within the function.
#' @param by_cluster Logical. If \code{TRUE}, density plots are faceted by individual and time-series variable.
#' If \code{FALSE}, only time-series variables are used for faceting. Default is \code{FALSE}.
#' @param cluster_ids Optional vector of cluster IDs to include in the plot. If \code{NULL}, all IDs are shown.
#' @param draw_person_pars Logical. If \code{TRUE}, samples are generated using person-specific parameters (random effects).
#' If \code{FALSE}, only population-level parameters are used. Defaults to \code{FALSE}.
#' @param n_draws Integer. Number of posterior draws to use for generating replicated datasets. Defaults to 20.
#' Ignored if \code{draws} is specified.
#' @param draws Optional vector of indices specifying which posterior draws to use. If \code{NULL}, \code{n_draws} samples are drawn randomly.
#' @param add_y_obs Logical. Whether to include the observed data distribution in the plot. Defaults to \code{TRUE}.
#' @param model_lab Optional character vector with labels for each model in \code{fit_list}. If \code{NULL}, defaults to "Model 1", "Model 2", etc.
#' @param y_rep_col Optional vector of colors for the posterior predictive densities of each model. If \code{NULL}, a default color palette is used.
#' @param y_obs_col Color for the observed data distribution. Default is \code{"#009E73"}.
#' @param y_rep_lw Line width of the observed data density curve. Default is \code{0.5}.
#' @param y_obs_lw Line width of the observed data density curve. Default is \code{1.1}.
#' @param y_rep_alpha Alpha transparency for the predictive density curves. Default is \code{0.5}.
#'
#' @details
#' This function performs graphical posterior predictive checks by overlaying kernel density estimates of replicated
#' data from the posterior with the observed data. This can be used to visually assess how well a fitted model captures
#' key distributional aspects of the observed time series. If \code{fit_list} is specified, multiple models can be
#' compared side-by-side in the same plot.
#'
#' If \code{draw_person_pars = TRUE}, simulated datasets incorporate subject-specific effects (random effects).
#' This requires that \code{monitor_person_pars = TRUE} was set during model fitting.
#'
#' @return A \code{ggplot} object showing density curves of observed and replicated data across time-series variables
#' (and optionally across individuals).
#'
#' @examples
#' \dontrun{
#' # Set up censored AR(1) model
#' ar1_cens <- mlts_model(q = 1, censor_left = -1)
#'
#' # Simulate data under the censored AR(1) model
#' simData <- mlts_sim(model = ar1_cens, N = 50, TP =100, default = TRUE)
#'
#' # Fit the model
#' fit_censAR <- mlts_fit(model = ar1_cens, data = simData$data,
#'                        id = "ID", ts = "Y1_cens", monitor_person_pars = TRUE)
#'
#' # As a comparison fit AR(1) model to the same data
#' ar1 <- mlts_model(q = 1)
#' fit_AR <- mlts_fit(model = ar1, data = simData$data,
#'                  id = "ID", ts = "Y1_cens", monitor_person_pars = T)
#'
#' # Run posterior predictive check
#' mlts_pp_check(fit_list = list(fit_censAR, fit_AR),
#'               model_lab = c("Censored AR(1)", "censored AR(1)"),
#'               y_rep_col = c("steelblue", "darkred"))
#'
#' # Check only selected individuals
#' mlts_pp_check(fit = fit_censAR, cluster_ids = c(1, 5, 10), by_cluster = TRUE)
#'
#' }
#'
#' @seealso \code{\link{mlts_posterior_sample}} for generating replicated data samples.
#' @export


mlts_pp_check <- function(
    fit,
    fit_list = NULL,
    ts = NULL,
    y_reps = NULL,
    by_cluster = FALSE,
    cluster_ids = NULL,
    draw_person_pars = FALSE,
    n_draws = 20,
    draws = NULL,
    add_y_obs = TRUE,
    model_lab = NULL,
    y_rep_col = NULL,
    y_obs_col = "#009E73",
    y_obs_lw = 1.1,
    y_rep_lw = 0.5,
    y_rep_alpha = 0.5){

  # base palette
  cbPalette <- c("#999999","#E69F00","darkblue","#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  # base settings
  if( is.null(y_rep_col) ){y_rep_col = cbPalette}
  if( is.null(model_lab)&!is.null(fit_list) ){model_lab = paste0("Model ", 1:length(fit_list), ": Y_rep")}
  if( is.null(model_lab)& is.null(fit_list) ){model_lab = paste0("Y_rep")}
  if( is.null(ts)       & is.null(fit_list) ) {ts = fit$standata$ts}
  if( is.null(ts)       &!is.null(fit_list) ) {ts = fit_list[[1]]$standata$ts}

  if(is.null(fit_list)) {

    if(is.null(y_reps)){
      y_reps <- mlts_posterior_sample(
        fit = fit,
        draw_person_pars = draw_person_pars,
        n_draws = n_draws,
        draws = draws
      )
    }

    # turn into data.frame
    y_reps <- do.call(rbind, y_reps)

    # bring into longer format by ts variables
    y_reps_long <- lapply(ts, function(x){
      cbind.data.frame(
        "ID" = y_reps$ID,
        "Y_class" =  model_lab,
        "Rep_no" = y_reps$Y_rep,
        "ts" = x,
        "Y" = y_reps[,x])
    })
    y_reps_long <- do.call(rbind, y_reps_long)


  } else {
    fit <- fit_list[[1]]

    y_reps_list <- list()
    for(j in 1:length(fit_list)){
      y_reps_list[[j]] <- mlts_posterior_sample(
        fit = fit_list[[j]],
        draw_person_pars = draw_person_pars,
        n_draws = n_draws,
        draws = draws
      )

      # turn into data.frame
      y_reps_list[[j]] <- do.call(rbind, y_reps_list[[j]])

      # bring into longer format by ts variables
      y_reps_list[[j]] <- lapply(ts, function(x){
        cbind.data.frame(
          "ID" = y_reps_list[[j]]$ID,
          "Y_class" = model_lab[j],
          "Rep_no" = paste0(model_lab[j], y_reps_list[[j]]$Y_rep),
          "ts" = x,
          "Y" = y_reps_list[[j]][,x])
      })

      y_reps_list[[j]] <- do.call(rbind, y_reps_list[[j]])

      }

    y_reps_long <- do.call(rbind, y_reps_list)

  }


  if(add_y_obs == TRUE){
    # bring observations into same format as replications
    y_obs_long <- lapply(ts, function(x){
      cbind.data.frame(
        "ID" = fit$data$num_id,
        "Y_class" = "Y_obs",
        "Rep_no" = 100000+1,
        "ts" = x,
        "Y" = fit$data[,x])
    })
    y_obs_long <- do.call(rbind, y_obs_long)

    y_reps_long <- rbind(y_obs_long, y_reps_long)
  }

  # data to plot
  if(!is.null(cluster_ids)){
    y_reps_long <- y_reps_long[y_reps_long$ID %in% cluster_ids,]
  }

  # set factor levels
  y_reps_long$Y_class = factor(y_reps_long$Y_class, levels = c("Y_obs", model_lab))

  # plot results
  pp.plot <- ggplot2::ggplot(y_reps_long, ggplot2::aes(x = .data$Y, group = .data$Rep_no, color = .data$Y_class)) +
    ggplot2::geom_density(alpha = y_rep_alpha,linewidth=y_rep_lw) +
    ggplot2::scale_color_manual(" ", values = c(y_obs_col,y_rep_col)) +
    ggplot2::theme_classic()

  if(add_y_obs==TRUE){
    pp.plot <- pp.plot +
      ggplot2::geom_density(
        inherit.aes = F, ggplot2::aes(x =  .data$Y, group = .data$Rep_no, color = .data$Y_class),
        color = y_obs_col,
        data=y_reps_long[y_reps_long$Y_class=="Y_obs",], alpha = 1, linewidth=y_obs_lw)
  }

  if(by_cluster == TRUE){
    pp.plot <- pp.plot +
      ggplot2::facet_wrap(ID~ts)
  } else {
    pp.plot <- pp.plot +
      ggplot2::facet_grid(~ts)
  }




  return(pp.plot)

}
