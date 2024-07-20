#' Plot results of mlts
#'
#' @param fit An object of class `mlts.fit`
#' @param type Type of plot.
#'             type = "fe" (Default)
#'                  Forest-plot of model coefficients.
#'             type = "re"
#'                  Plot of individual (random) effects
#'             type = "int"
#'                  Experimental: Plot within-level interactions.
#'             type = "re.cor"
#'                  Combined plot depicting the distribution of individual parameter
#'                  estimates (posterior summary statistics as provided by `bpe`), as
#'                  well as bivariate scatter plots.
#' @param what Character. For `type = "fe"`, indicate which parameters should be
#' included in the plot by setting `what` to "all" (the default), or one (or multiple)
#' of "Fixed effect", "Random effect SD", "RE correlation", "Outcome prediction", "RE prediction",
#' "Item intercepts", "Loading", or "Measurement Error SD".
#' @param bpe The Bayesian point estimate is, by default, the median of the
#' posterior distribution (`bpe = "median"`). Set `bpe = "mean"` to use
#' the mean of the posterior distribution as point estimates.
#' @param sort_est Add parameter label for sorting of random effects.
#' @param xlab Title for the x axis.
#' @param ylab Title for the y axis.
#' @param facet_ncol Number of facet columns (see `ggplot2::facet_grid`).
#' @param dot_size numeric, size of the dots that indicate the point estimates.
#' @param dot_color character. indicating the color of the point estimates.
#' @param dot_shape numeric. shape of the dots that indicate the point estimates.
#' @param errorbar_color character. Color of error bars.
#' @param errorbar_width integer. Width of error bars.
#' @param add_true logical. If model was fitted with simulated data using `mlts_sim`,
#' true population parameter values can be plotted as reference by setting the argument ot `TRUE`.
#' @param true_color character. Color of points depicting true population parameter used in the data generation.
#' @param true_shape integer. Shape of points depicting true population parameter used in the data generation.
#' @param true_size integer. Size of points depicting true population parameter used in the data generation.
#' @param hide_xaxis_text logical. Hide x-axis text if set to `TRUE`.
#' @param par_labels character vector. User-specified labels for random effect parameters
#' can be specified.
#' @param labels_as_expressions logical. Should parameter names on plot labels be printed
#' as mathematical expressions? Defaults to `FALSE`. Still experimental.
#' @return Returns a `ggplot`-object .
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
#'   time = "time",
#'   tinterval = 1 # interval for approximation of continuous-time dynamic model,
#' )
#'
#' # inspect model summary
#' mlts_plot(fit, type = "fe", what = "Fixed effect")
#' }
mlts_plot <- function(fit, type = c("fe", "re", "re.cor", "int"), bpe = c("median", "mean"),
                      what = c("all", "Fixed effect", "Random effect SD", "RE correlation",
                               "Outcome prediction", "RE prediction", "Item intercepts",
                               "Loading", "Measurement Error SD"),
                      sort_est = NULL, xlab = NULL, ylab = NULL,
                      facet_ncol = 1, dot_size = 1, dot_color = "black", dot_shape = 1,
                      errorbar_color = "black", errorbar_width = 0.3,
                      add_true = FALSE, true_color = "red", true_shape = 22, true_size = 1,
                      hide_xaxis_text = TRUE,
                      par_labels = NULL, labels_as_expressions = FALSE){

  what <- match.arg(what, several.ok = TRUE)
  type <- match.arg(type)
  bpe <- match.arg(bpe)



  if(type == "fe"){ # fixed effect estimates
    if(bpe == "median"){
      fit$pop.pars.summary$bpe = fit$pop.pars.summary$`50%`
      xlab <- ifelse(!is.null(xlab), xlab, "Posterior Median (95% Credibility Interval)")
    } else {
      fit$pop.pars.summary$bpe = fit$pop.pars.summary$mean
      xlab <- ifelse(!is.null(xlab), xlab, "Posterior Mean (95% Credibility Interval)")
    }

    # extract data used for plotting
    p.data <- fit$pop.pars.summary
    # order fix effects
    p.data$Param <- factor(p.data$Param, levels = rev(p.data$Param))
    # add parameter type used for facets
    p.data$Param_type <- p.data$Type

    # subset data by what
    if(!("all" %in% what)){
      p.data <- p.data[p.data$Type %in% c(what),]
    }

    # order parameter type
    p.data$Param_type <- factor(p.data$Param_type, levels = c(
      "Fixed effect", "Random effect SD", "RE correlation",
      "Outcome prediction", "RE prediction",
      "Item intercepts", "Loading", "Measurement Error SD"
    ))



    # build general plot
    aes <- ggplot2::aes
    P <- ggplot2::ggplot(data = p.data, aes(y = .data$Param, x = .data$bpe)) +
          ggplot2::geom_point(color = dot_color, size = dot_size, shape = dot_shape) +
          ggplot2::geom_errorbar(aes(xmin = .data$`2.5%`, xmax = .data$`97.5%`),
                                 color = errorbar_color,
                                 width = errorbar_width) +
          ggplot2::facet_wrap(~.data$Param_type, ncol = facet_ncol, scales = "free", shrink = TRUE) +
          ggplot2::labs(x = xlab, y = ylab) +
          ggplot2::scale_x_continuous(n.breaks = 8) +
          ggplot2::theme_bw()

    # optional add true scores
    if(add_true == TRUE){
    P <- P + ggplot2::geom_point(aes(x = .data$true.val),
                                 color = true_color,
                                 shape = true_shape,
                                 size = true_size)
    }
  }

  if(type == "re"){ # random effect estimates
    if(is.null(nrow(fit$person.pars.summary))){
      stop("No individual (random) parameter results available. Model needs to be refitted with monitor.person.pars = TRUE in VARfit().")
    }
    xlab <- ifelse(!is.null(xlab), xlab, colnames(fit$person.pars.summary)[2])

    if(bpe == "median"){
      fit$person.pars.summary$bpe = fit$person.pars.summary$`50%`
      ylab <- ifelse(!is.null(ylab), ylab, "Posterior Median (95% Credibility Interval)")
    } else {
      fit$person.pars.summary$bpe = fit$person.pars.summary$mean
      ylab <- ifelse(!is.null(ylab), ylab, "Posterior Mean (95% Credibility Interval)")
    }

    # evaluate model
    infos <- mlts_model_eval(fit$model)

    # extract data used for plotting
    p.data <- fit$person.pars.summary

    # store original column names as a separate variable
    p.data$ID <- p.data[,2]

    # if sorting is requested, change levels of ID variable
    if(!is.null(sort_est)){
      sorted = p.data[p.data$Param == sort_est,]
      sorted = sorted[order(sorted$mean),]
      p.data$ID <- factor(p.data$ID, levels = sorted$ID)
    } else {
      p.data$ID = factor(p.data$ID, levels = sort(unique(p.data$ID)))
    }


    # order random effects
    p.data$Param = factor(p.data$Param, levels = infos$re_pars$Param)

    # build general plot
    aes = ggplot2::aes
    P <- ggplot2::ggplot(data = p.data, ggplot2::aes(x = .data$ID, y = .data$mean)) +
           ggplot2::geom_point() +
           ggplot2::geom_errorbar(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`),
                                  color = errorbar_color,
                                  width = errorbar_width) +
           ggplot2::facet_wrap(~.data$Param, scale = "free", ncol = facet_ncol) +
           ggplot2::theme_bw() +
           ggplot2::labs(y = ylab, x = xlab) +
           ggplot2::theme(axis.text.x = ggplot2::element_text(
             angle = 45, hjust = 1, vjust = 1, size = 7))

     # optional add true scores
     if(add_true == TRUE){
     P <- P + ggplot2::geom_point(aes(y = .data$true.val),
                                  color = true_color,
                                  shape = true_shape,
                                  size = true_size)
     }

     if(hide_xaxis_text){
     P <- P + ggplot2::theme(axis.text.x = ggplot2::element_blank())
     }

  }

  if(type == "re.cor"){ # random effect estimates

    # create a mixture plot depicting the distribution of random effect estimates
    # as histogram on the plot-grid diagnonal as well as scatterplots on the off-diagonal
    # panels

    # evaluate model
    infos <- mlts_model_eval(fit$model)

    # get the data
    if(bpe == "median"){
      fit$person.pars.summary$bpe = fit$person.pars.summary$`50%`
      ylab <- ifelse(!is.null(ylab), ylab, "Posterior Median (95% Credibility Interval)")
    } else {
      fit$person.pars.summary$bpe = fit$person.pars.summary$mean
      ylab <- ifelse(!is.null(ylab), ylab, "Posterior Mean (95% Credibility Interval)")
    }

    # get the labels of included RE params
    p.data = fit$person.pars.summary

    if(!is.null(par_labels)){
      p.data$Param = factor(p.data$Param, levels = par_labels,
                            labels = names(par_labels))
    } else {
      p.data$Param = factor(p.data$Param,
                            levels = infos$re_pars$Param,
                            labels = infos$re_pars$Param)
    }

    re_pars = levels(p.data$Param)
    n_re_pars = length(re_pars)
    # loop plotting over all combinations of re_pars
    re_par_combi = data.frame(
      "RE1" = rep(re_pars, each = n_re_pars),
      "Re2" = rep(re_pars, times = n_re_pars)
    )
    p_list = list()
    sub = c()

    suppressWarnings({
      suppressMessages({

    for(i in 1:nrow(re_par_combi)){

      sub = p.data[p.data$Param %in% unlist(re_par_combi[i,]),]

      if(length(unique(sub$Param)) == 1){
        sub$lab_x = re_par_combi[i,1]
        sub$lab_y = re_par_combi[i,1]
        p_list[[i]] <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$bpe)) +
          ggplot2::geom_histogram(fill = "grey", color = "black") +
          ggplot2::theme_classic() +
          ggplot2::scale_x_continuous(n.breaks = 7) +
          ggplot2::theme(strip.placement = "outside",
                strip.background = ggplot2::element_blank(),
                strip.text = ggplot2::element_text(face = 2),
                axis.title = ggplot2::element_blank())

        if(labels_as_expressions == T){
          p_list[[i]] <- p_list[[i]] +
            ggplot2::facet_grid(cols = ggplot2::vars(.data$lab_x), switch = "y",
                       labeller = ggplot2::label_parsed)
        } else {
          p_list[[i]] <- p_list[[i]] +
            ggplot2::facet_grid(cols = ggplot2::vars(.data$lab_x), switch = "y")
        }

      } else {
        df = data.frame(
          y = sub$bpe[sub$Param == re_par_combi[i,1]],
          x = sub$bpe[sub$Param == re_par_combi[i,2]],
          lab_x = re_par_combi[i,2],
          lab_y = re_par_combi[i,1]
        )
        p_list[[i]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
          ggplot2::geom_point(color = dot_color, shape = dot_shape, fill = "grey") +
          ggplot2::theme_classic() +
          ggplot2::stat_smooth(method = "lm", se = FALSE) +
          ggplot2::scale_x_continuous(n.breaks = 7) +
          ggplot2::scale_y_continuous(n.breaks = 7) +
          ggplot2::theme(strip.placement = "outside",
                strip.background = ggplot2::element_blank(),
                strip.text = ggplot2::element_text(face = 2),
                axis.title = ggplot2::element_blank()
          )

  #       if(labels_as_expressions == T){
  #         p_list[[i]] <- p_list[[i]] +
  #           ggplot2::facet_grid(#cols = ggplot2::vars(lab_x),
  # #                              rows = ggplot2::vars(lab_y),
  #                               switch = "y",
  #                      labeller = ggplot2::label_parsed)
  #       } else {
  #         p_list[[i]] <- p_list[[i]] +
  #           ggplot2::facet_grid(#cols = ggplot2::vars(lab_x),
  #                             #  rows = ggplot2::vars(lab_y),
  #                               switch = "y")
  #       }
      }
    }

    P = cowplot::plot_grid(plotlist = p_list, align = "hv")

      })
    })

  }


  if(type == "int"){

    stop("Work in progress!")

    # check if interactions are involved
    if(fit$standata$n_int==0){
      stop("Model does not contain interaction effects on the dynamic within level.")
    }

    # evaluate model
    infos <- mlts_model_eval(fit$model)

    # Interaction plot
    ## unstandardized effects

    ## 1. calculate predicted scores of y using the range of the observed scores
    int_pars = infos$fix_pars_dyn[infos$fix_pars_dyn$isINT == 1,]
    int_y = as.integer(int_pars$Dout)
    int_x1 = as.integer(int_pars$Dpred)
    int_x2 = as.integer(int_pars$Dpred2)

    # range of
#    x_axis =

    ## for now: only plot standardized effects
#    mlts_standardized(fit, what = "both")
    x_axis = seq(-3,3, by = 0.01)


    # general setup for interaction plot
    # 1. extract samples
#    b_samples =

    # 2. specify range to plot
    x_axis <- seq(-3,3, by = 0.01)

    # 3. calculate predicted scores
    # pred_mean <-  b_std[,6] %*% t(x_axis)
    # pred_SDup <-  b_std[,6] %*% t(x_axis) + b_std[,7]%*% t(x_axis)
    # pred_SDdown <-  b_std[,6] %*% t(x_axis) - b_std[,7]%*% t(x_axis)
    #
    # fit_mean <- as.vector(colMeans(pred_mean))
    # fit_SDup <- as.vector(colMeans(pred_SDup))
    # fit_SDdown <- as.vector(colMeans(pred_SDdown))
    #
    # tmp <-data.frame(rep(x_axis,3),c(fit_SDup,fit_mean,fit_SDdown),c(rep("+1SD",length(x_axis)),rep("mean",length(x_axis)),rep("-1SD",length(x_axis))) )
    # names(tmp) <- c("Mindful_attention","Negative_Affect","PA_previous")
    #
    # fit_mean_2.5 <- as.vector(apply(pred_mean,2,function(x) quantile(x, probs=0.025)))
    # fit_mean_97.5 <- as.vector(apply(pred_mean,2,function(x) quantile(x, probs=0.975)))
    # fit_SDupT_2.5 <- as.vector(apply(pred_SDup,2,function(x) quantile(x, probs=0.025)))
    # fit_SDupT_97.5 <- as.vector(apply(pred_SDup,2,function(x) quantile(x, probs=0.975)))
    # fit_SDdownT_2.5 <- as.vector(apply(pred_SDdown,2,function(x) quantile(x, probs=0.025)))
    # fit_SDdownT_97.5 <- as.vector(apply(pred_SDdown,2,function(x) quantile(x, probs=0.975)))
    # tmp$CI2.5 <- c(fit_SDupT_2.5,fit_mean_2.5,fit_SDdownT_2.5)
    # tmp$CI97.5 <- c(fit_SDupT_97.5,fit_mean_97.5,fit_SDdownT_97.5)






  }



  return(P)

}
