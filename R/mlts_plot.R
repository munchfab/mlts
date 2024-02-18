#' Plot results of mlts
#'
#' @param fit An object of class `mlts.fit`
#' @param type Type of plot.
#'             type = "fixef" (Default)
#'                  Forest-plot of model coefficients.
#'             type = "ranef"
#'                  Plot of individual (random) effects
#'             type = "raneftab"
#'                  (not yet implemented)
#' @param bpe The Bayesian point estimate is, by default, the median of the
#' posterior distribution (`bpe = "median"`). Set `bpe = "mean"` to use
#' the mean of the posterior distribution as point estimates.
#' @param sort.est Add parameter label for sorting of random effects.
#' @param xlab Title for the x axis.
#' @param ylab Title for the y axis.
#' @param facet.ncol Number of facet columns (see `ggplot2::facet_grid`).
#' @param dot.size Numeric, size of the dots that indicate the point estimates.
#' @param dot.color Character vector, indicating the color of the point estimates.
#' @param dot.shape Numeric, shape of the dots that indicate the point estimates.
#'
#' @return Returns a `ggplot`-object .
#'
#' @export
#'
#' @examples TBA
mlts_plot <- function(fit, type = c("fe", "re", "re.cor"), bpe = c("median", "mean"),
                      sort.est = NULL, xlab = NULL, ylab = NULL,
                      facet.ncol = 1, dot.size = 1, dot.color = "black", dot.shape = 1,
                      err.bar.col = "black", err.bar.width = 0.3,
                      add.true = FALSE, true.col = "red", true.shape = 22, true.size = 1,
                      hide.xaxis.text = T,
                      par_labels = NULL, labels.as.expression = F){

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
    p.data$Param_type <- fit$VARmodel$Type
    # order parameter type
    p.data$Param_type <- factor(p.data$Param_type, levels = c(
      "Fix effect", "Random effect SD", "RE correlation",
      "Outcome prediction", "RE prediction",
      "Item intercepts", "Loading", "Measruement Error SD"
    ))


    # build general plot
    aes <- ggplot2::aes
    P <- ggplot2::ggplot(data = p.data, aes(y = Param, x = bpe)) +
          ggplot2::geom_point(color = dot.color, size = dot.size, shape = dot.shape) +
          ggplot2::geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`),
                                 color = err.bar.col,
                                 width = err.bar.width) +
          ggplot2::facet_wrap(~Param_type, ncol = facet.ncol, scales = "free", shrink = T) +
          ggplot2::labs(x = xlab, y = ylab) +
          ggplot2::scale_x_continuous(n.breaks = 8) +
          ggplot2::theme_bw()

    # optional add true scores
    if(add.true == TRUE){
    P <- P + ggplot2::geom_point(aes(x = true.val),
                                 color = true.col,
                                 shape = true.shape,
                                 size = true.size)
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
    infos <- VARmodelEval(fit$VARmodel)

    # extract data used for plotting
    p.data <- fit$person.pars.summary

    # store original column names as a separate variable
    p.data$ID <- p.data[,2]

    # if sorting is requested, change levels of ID variable
    if(!is.null(sort.est)){
      sorted = p.data[p.data$Param == sort.est,]
      sorted = sorted[order(sorted$mean),]
      p.data$ID <- factor(p.data$ID, levels = sorted$ID)
    } else {
      p.data$ID = factor(p.data$ID, levels = sort(unique(p.data$ID)))
    }


    # order random effects
    p.data$Param = factor(p.data$Param, levels = infos$re_pars$Param)

    # build general plot
    aes = ggplot2::aes
    P <- ggplot2::ggplot(data = p.data, aes(x = ID, y = mean)) +
           ggplot2::geom_point() +
           ggplot2::geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),
                                  color = err.bar.col,
                                  width = err.bar.width) +
           ggplot2::facet_wrap(~Param, scale = "free", ncol = facet.ncol) +
           ggplot2::theme_bw() +
           ggplot2::labs(y = ylab, x = xlab) +
           ggplot2::theme(axis.text.x = ggplot2::element_text(
             angle = 45, hjust = 1, vjust = 1, size = 7))

     # optional add true scores
     if(add.true == TRUE){
     P <- P + ggplot2::geom_point(aes(y = true.val),
                                  color = true.col,
                                  shape = true.shape,
                                  size = true.size)
     }

     if(hide.xaxis.text){
     P <- P + ggplot2::theme(axis.text.x = ggplot2::element_blank())
     }

  }

  if(type == "re.cor"){ # random effect estimates

    # create a mixture plot depicting the distribution of random effect estimates
    # as histogram on the plot-grid diagnonal as well as scatterplots on the off-diagonal
    # panels

    # evaluate model
    infos <- VARmodelEval(fit$VARmodel)

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
    for(i in 1:nrow(re_par_combi)){

      sub = p.data[p.data$Param %in% unlist(re_par_combi[i,]),]

      if(length(unique(sub$Param)) == 1){
        sub$lab_x = re_par_combi[i,1]
        sub$lab_y = re_par_combi[i,1]
        p_list[[i]] <- ggplot(sub, aes(x = bpe)) +
          geom_histogram(fill = "grey", color = "black") +
          theme_classic() +
          scale_x_continuous(n.breaks = 7) +
          theme(strip.placement = "outside",
                strip.background = element_blank(),
                strip.text = element_text(face = 2),
                axis.title = element_blank())

        if(labels.as.expression == T){
          p_list[[i]] <- p_list[[i]] +
            facet_grid(cols = vars(lab_x), rows = vars(lab_y), switch = "y",
                       labeller = label_parsed)
        } else {
          p_list[[i]] <- p_list[[i]] +
            facet_grid(cols = vars(lab_x), rows = vars(lab_y), switch = "y")
        }

      } else {
        df = data.frame(
          y = sub$bpe[sub$Param == re_par_combi[i,1]],
          x = sub$bpe[sub$Param == re_par_combi[i,2]],
          lab_x = re_par_combi[i,2],
          lab_y = re_par_combi[i,1]
        )
        p_list[[i]] <- ggplot(df, aes(x = x, y = y)) +
          geom_point(color = dot.color, shape = dot.shape, fill = "grey") +
          theme_classic() +
          stat_smooth(method = "lm", se = F) +
          scale_x_continuous(n.breaks = 7) +
          scale_y_continuous(n.breaks = 7) +
          theme(strip.placement = "outside",
                strip.background = element_blank(),
                strip.text = element_text(face = 2),
                axis.title = element_blank()
          )

        if(labels.as.expression == T){
          p_list[[i]] <- p_list[[i]] +
            facet_grid(cols = vars(lab_x), rows = vars(lab_y), switch = "y",
                       labeller = label_parsed)
        } else {
          p_list[[i]] <- p_list[[i]] +
            facet_grid(cols = vars(lab_x), rows = vars(lab_y), switch = "y")
        }
      }
    }

    P = cowplot::plot_grid(plotlist = p_list, align = "hv")

  }


  return(P)

}
