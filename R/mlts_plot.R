#' Plot results of mlts
#'
#' @param fit An object of class `mlts.fit`
#' @param what one of "fixef", "ranef", "raneftab"
#' @return A `list` that can be passed to `stan()`.
#' @export
#'
#' @examples TBA
mlts_plot <- function(fit, type = c("fe", "re", "re.cor"), sort.est = NULL, xlab = NULL, ylab = NULL,
                      n.facet.col = 1,pt.size = 1, pt.col = "black", pt.shape = 1,
                      err.bar.col = "black", err.bar.width = 0.3,
                      add.true = FALSE, true.col = "red", true.shape = 22, true.size = 1){

  type <- match.arg(type)

  if(type == "fe"){ # fixed effect estimates
    xlab <- ifelse(!is.null(xlab), xlab, "Posterior Mean (95% Credibility Interval)")
    ylab <- ifelse(!is.null(ylab), ylab, "Parameter")

    # extract data used for plotting
    p.data <- fit$pop.pars.summary
    # order fix effects
    p.data$Param <- factor(p.data$Param, levels = rev(p.data$Param))
    # add parameter type used for facets
    p.data$Param_type <- VARresults$VARmodel$Type

    # build general plot
    aes <- ggplot2::aes
    P <- ggplot2::ggplot(data = p.data, aes(y = Param, x = mean)) +
          ggplot2::geom_point(color = pt.col, size = pt.size, shape = pt.shape) +
          ggplot2::geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`),
                                 color = err.bar.col,
                                 width = err.bar.width) +
          ggplot2::facet_wrap(~Param_type, ncol = n.facet.col, scales = "free") +
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
    xlab <- ifelse(!is.null(xlab), xlab, colnames(fit$person.pars.summary)[2])
    ylab <- ifelse(!is.null(ylab), ylab, "Posterior Mean (95% Credibility Interval)")

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
           ggplot2::facet_wrap(~Param, scale = "free", ncol = n.facet.col) +
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
  }

  if(type == "re.cor"){ # random effect estimates


    ###### ADD LATER :::::::::::::::::

    # # get person par data
    # data = fit$person.pars.summary
    #
    # # prepare empty data frame
    # n_ids = max(fit$person.pars.summary$num_id)
    # re.pars = unique(fit$person.pars.summary$Param)
    # n_re.pars = length(re.pars)
    #
    #
    #
    # p.data = data.frame(
    #   "ID" = rep(1:n_ids, times = n_ids*n_re.pars*n_re.pars),
    #   "Param_x" = rep(unique(fit$person.pars.summary$Param), each = n_ids,
    #                   times = n_re.pars),
    #   "Param_y" = rep(unique(fit$person.pars.summary$Param), each = n_ids*n_re.pars),
    #   "x" = NA,
    #   "y" = NA
    # )
    #
    # # add values
    # for(i in 1:nrow(p.data)){
    #   p.data$x[i] = data$mean[data$num_id==p.data$ID[i] & data$Param == p.data$Param_x[i]]
    #   p.data$y[i] = data$mean[data$num_id==p.data$ID[i] & data$Param == p.data$Param_y[i]]
    # }
    #
    # # plot
    # ggplot(p.data, aes(x = x, y = y)) +
    #   geom_point() +
    #   facet_grid(Param_x ~ Param_y, scales = "free_x")
  }


  return(P)

}
