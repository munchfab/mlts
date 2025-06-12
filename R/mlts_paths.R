#' Plot Paths for Two-Level VAR Model
#'
#' @description
#' The `mlts_paths` function depcits models specified using `mlts_model` as a
#' path diagram.
#'
#' @param model `data.frame`. Output of \code{\link[mlts]{mlts_model}} and related functions.
#' @param asp_decomp A numeric value specifying the aspect ratio for the decomposition
#'   plot region. Defaults to 0.25.
#' @param asp_w_b A numeric value specifying the aspect ratio between the within-level
#'   and between-level sections. Defaults to 0.5.
#' @param fig_margins.x A numeric vector of length 2 defining the horizontal margins
#'   of the plot. Defaults to `c(0, 8)`.
#' @param fig_margins.y A numeric vector of length 2 defining the vertical margins
#'   of the plot. Defaults to `c(0, 8)`.
#' @param width Width of the plot in inches. Defaults to 7.
#' @param height Height of the plot in inches. Defaults to 6.
#' @param file A character string specifying the path to save the plot. If `NULL`,
#'   the plot will not be saved. Defaults to `NULL`.
#' @param asp The overall aspect ratio of the plot, computed as `height / width`.
#'   Defaults to `height / width`.
#' @param family Font family used in the plot. Defaults to `"serif"`.
#' @param cex_b Numeric value specifying the scaling of text in the between-level
#'   section. Defaults to 0.8.
#' @param cex_w Numeric value specifying the scaling of text in the within-level
#'   section. Defaults to 0.8.
#' @param cex_decomp Numeric value specifying the scaling of text in the decomposition
#'   section. Defaults to 0.8.
#' @param cex_loads Numeric value specifying the scaling of text of loading parameters.
#'   Defaults to 0.8.
#' @param adj_load_x Numeric value specifying the x-axis offset loading parameter labels.
#'   Defaults to 1.25.
#' @param y_ind_labs A vector of character strings with names of observed variables.
#' @param y_fac_labs A vector of character strings with factor labels to replace numeric indices in parameter names.
#' @param y_fac_lab_sep A character string to separate multiple factor labels. Defaults to ",".
#' @param remove_lag_lab Logical. Remove lag index from phi-parameter labels. Defaults to `FALSE`.
#' @param b_style A character string specifying the style of the between-level plot
#'   ("h" for horizontal). Defaults to `"h"`.
#' @param arrHead_w Numeric values controlling the arrowhead size for
#'   within-level paths. Defaults to 0.16.
#' @param arrHead_b Numeric values controlling the arrowhead size for
#'   between-level paths. Defaults to 0.16.
#' @param scale_decomp_ind Numeric. Specify the scaling factor for manifest indicators in the
#'   decomposition section.
#' @param scale_decomp_F Numeric. Specify the scaling factor for latent factors in the
#'   decomposition section.
#' @param scale_within Numeric. Specify the scaling factor for latent factors in the
#'   within-level section.
#' @param scale_within_inno Numeric. Specify the scaling factor for innovations in the
#'   within-level section.
#' @param scale_between Numeric. Specify the scaling factor for factors in the
#'   between-level section.
#' @param scale_between Numeric. Specify the scaling factor for factors in the
#'   between-level section.
#' @param scale_int Numeric. Specify the scaling factor for interaction factors in the
#'   within-level section.
#' @param lwd_nodes Line width for node borders in the plot. Defaults to 1.7.
#' @param rand_dot_pos Numeric value controlling the random dot position in the plot.
#'   Defaults to 0.5.
#' @param res The nominal resolution in ppi. Defaults to 320.
#' @param units A character string specifying the units for saving the plot. Defaults to "in".
#' @param pointsize Numeric value specifying the font point size for the plot. Defaults to 10.
#' @param type A character string specifying the file type for the saved plot (e.g., "cairo"). Defaults to "cairo".
#' @param ... Additional arguments passed to internal plotting functions.
#'
#' @details
#' This function calculates positions, radii, and labels for nodes and arrows based
#' on the model structure and its parameters. It divides the plot into sections:
#'
#' - **Decomposition**: Shows the breakdown of observed variables into within- and
#'   between-level components.
#' - **Within-Level Dynamics**: Illustrates autoregressive and cross-lagged paths
#'   between variables at the within level.
#' - **Between-Level Dynamics**: Depicts random effects, covariates, and their
#'   interrelations at the between level.
#'
#' Depending on the model structure (e.g., maximum lag, number of random effects,
#' presence of interaction terms), the function dynamically adjusts the visualization.
#'
#' @return A graphical object representing the path diagram of the model.
#'
#' @examples
#' # A two-level second-order autoregressive model
#' model <- mlts_model(q = 1, max_lag = 2)
#'
#' # Plot the paths
#' mlts_paths(model)
#'
#' @export

mlts_paths <- function(
    model,
    asp_decomp = 0.25,
    asp_w_b = 0.5,
    fig_margins.x = c(0,8),
    fig_margins.y = c(0,8),
    width = 7,
    height = 6,
    file = NULL,
    asp = height/width,
    family = "serif",
    cex_b = 0.8,
    cex_w = 1,
    cex_decomp = 1,
    cex_loads = 0.8,
    b_style = "h",
    arrHead_w = 0.16,
    arrHead_b = 0.16,
    scale_decomp_ind = 0.35,
    scale_decomp_F = 0.45,
    scale_within = 0.3,
    scale_within_inno = 0.2,
    scale_between = 0.3,
    scale_int = 0.25,
    lwd_nodes = 1.7,
    rand_dot_pos = 0.4,
    units = "in",
    res = 700,
    pointsize = 10,
    type = "cairo",
    # new options
    y_ind_labs = NULL,
    y_fac_labs = NULL, # + cex_loads
    y_fac_lab_sep = ",",
    remove_lag_lab = FALSE,
    adj_load_x = 1.25,
    ...
    ) {

# get infos on model
infos <- mlts_model_eval(model)

# set the widths between nodes in the within-level model:
if(infos$maxLag == 1){
  within_x_pos <- c(60,15)
  within_inno_pos <- 80
}
if(infos$maxLag == 2){
  within_x_pos <- c(65,40,15)
  within_inno_pos <- 80
}

# adjust some inputs to create sensible default plots:
if(infos$n_random > 6){
  scale_between = scale_between + 0.1
}
if(infos$q == 1){
  scale_within = scale_within + 0.1
}


# decomposition ================================================================
begin.decomp.x <- fig_margins.x[1]
end.decomp.x <- fig_margins.x[1] + diff(fig_margins.x) * asp_decomp
w.decomp.x <- end.decomp.x - begin.decomp.x
w.decomp.y <- fig_margins.y[2] - fig_margins.y[1]
# within =======================================================================
begin.wth.x <- fig_margins.x[1] + end.decomp.x
end.wth.x <- fig_margins.x[2]
w.wth.x <- end.wth.x - begin.wth.x
begin.wth.y <- fig_margins.y[1] + asp_w_b*diff(fig_margins.y)
end.wth.y <- fig_margins.y[2]
w.wth.y <- end.wth.y - begin.wth.y
# between ======================================================================
begin.btw.x <- fig_margins.x[1] + end.decomp.x
end.btw.x <- fig_margins.x[2]
w.btw.x <- end.wth.x - begin.wth.x
begin.btw.y <- fig_margins.y[1]
end.btw.y <- fig_margins.y[1] + asp_w_b*diff(fig_margins.y)
w.btw.y <- end.btw.y - begin.btw.y




# DECOMPOSITION ================================================================
ind_nodes = mlts_path_decomp(
  infos = infos,
  asp = asp,
  begin.decomp.x = begin.decomp.x,
  end.decomp.x = end.decomp.x,
  w.decomp.y = w.decomp.y,
  fig_margins.y = fig_margins.y,
  scale_decomp_ind = scale_decomp_ind,
  scale_decomp_F = scale_decomp_F,
  y_ind_labs = y_ind_labs,
  y_fac_labs = y_fac_labs
)


# WITHIN =======================================================================
n_nod_w_y <- infos$q
n_nod_w <- infos$q
n_TP <- infos$maxLag:0
#w_nodes <- infos$fix_pars_dyn
# n_pos_w < n_TP
# calculate positions
w_poses_y = get_mid_points(n = n_nod_w, lims = c(end.wth.y, begin.wth.y))
#w_poses_x = get_mid_points(n = n_TP, lims = c(end.wth.x, begin.wth.x))
# get the nodes
w_nodes = data.frame()
for(i in 1:n_nod_w_y){
  for(j in n_TP){
    row = cbind.data.frame(
      "lag" = j,
      "time" = ifelse(j ==0, "italic(t)", paste0("(italic(t-",j,"))")),
      "construct" = i,
      "prefix" = ifelse(infos$isLatent == 1, "eta", "Y"),
      "w_cen" = 1
    )
    w_nodes = rbind(w_nodes, row)
  }
}
for(i in 1:nrow(w_nodes)){
  if(infos$is_wcen[w_nodes$construct[i]]==0){
    w_nodes$w_cen[i] = 0
  }
}

# create the labels
w_nodes$par <- ifelse(w_nodes$w_cen == 1,
        paste0(w_nodes$prefix,"[italic(", w_nodes$construct,"*i)*", w_nodes$time, "]^W"),
        paste0("Y[italic(", w_nodes$construct,"*i)*", w_nodes$time, "]"))
w_nodes$lab <- plotmath_labeller(x = w_nodes$par, y_fac_labs = y_fac_labs)
# get positions
w_ys <- get_mid_points(n = n_nod_w_y, lims = c(end.wth.y, begin.wth.y))
w_xs <- get_mid_points(n = 100, lims = c(begin.wth.x, end.wth.x))[within_x_pos]
w_radx <- diff(get_mid_points(n = 100, lims = c(begin.wth.x, end.wth.x))[c(10,30)])*scale_within
w_inno_x <- get_mid_points(n = 100, lims = c(begin.wth.x, end.wth.x))[within_inno_pos]
w_inno_x2 <- w_inno_x + w_radx
for(i in 1:nrow(w_nodes)){
  w_nodes$midx[i] <- w_xs[w_nodes$lag[i]+1]
  w_nodes$midy[i] <- w_ys[w_nodes$construct[i]]
}
# INNOVATION VARIANCES =========================================================
inno <- infos$fix_pars
inno <- inno[grepl(inno$Param_Label, pattern = "Innovation Variance"),]
inno$midx <- w_inno_x
inno$midy <- w_ys[infos$is_wcen == 1]
inno_radx <- diff(get_mid_points(n = 100, lims = c(begin.wth.x, end.wth.x))[c(10,30)])*scale_within_inno

for(i in 1:nrow(inno)){
  inno$q[i] <- as.integer(strsplit(inno$Param[i], split = "_")[[1]][2])
  if(!is.null(y_fac_labs)){
    inno$lab[i] <- paste0("zeta[",y_fac_labs[i],"*italic(it)]")
  } else {
    inno$lab[i] <- paste0("zeta[",inno$q[i],"*italic(it)]")
  }
  inno$x0[i] <- inno$midx[i] - inno_radx
  inno$x1[i] <- w_nodes$midx[w_nodes$construct==inno$q[i]&w_nodes$lag==0] + w_radx
  inno$y0[i] <- w_ys[inno$q[i]]
  inno$y1[i] <- w_ys[inno$q[i]]
}

# innovation covariance factors
if(infos$n_inno_covs == 1){
  inno_cov <- infos$fix_pars
  inno_cov <- inno_cov[grepl(inno_cov$Param_Label, pattern = "Covariance"),]
  inno_cov$midx <- w_inno_x
  inno_cov$midy <- get_mid_points(n = 1, w_ys[1:2])
  inno_radx <- diff(get_mid_points(n = 100, lims = c(begin.wth.x, end.wth.x))[c(10,30)])*scale_within_inno

  inno_cov = rbind(inno_cov, inno_cov)
  for(i in 1:2){
#    inno$q[i] <- as.integer(strsplit(inno$Param[i], split = "_")[[1]][2])
    inno_cov$lab <- paste0("eta[zeta*",12,"*italic(it)]")
    inno_cov$x0[i] <- inno_cov$midx[i] - inno_radx
    inno_cov$x1[i] <- w_nodes$midx[w_nodes$construct==inno$q[i]&w_nodes$lag==0]
    inno_cov$y0[i] <- inno_cov$midy[i]
    inno_cov$y1[i] <- w_ys[inno$q[i]]
    inno_cov$arr_pos[i] <- get_arr_pos(inno_cov$x0[i],inno_cov$x1[i],y0 = inno_cov$y0[i],inno_cov$y1[i],radx = inno_radx+0.5*arrHead_w)
    inno_cov$loads[i] <- infos$inno_cov_load[i]
    }
  }

# innovation correlations
if(nrow(infos$inno_cors)>0){
  inno_cor = infos$inno_cors
  # get involved innovations
  inno_cor$inno1 = substr(inno_cor$Param, nchar(inno_cor$Param)-1, nchar(inno_cor$Param)-1)
  inno_cor$inno2 = substr(inno_cor$Param, nchar(inno_cor$Param), nchar(inno_cor$Param))
  inno_cor$diff = abs(as.numeric(inno_cor$inno1)-as.numeric(inno_cor$inno2))

  # set curvature
  inno_cor$curve = ifelse(inno_cor$diff>2, 0.5, inno_cor$diff * 0.22)

  # get midpoints
  inno_cor$inno1_midy = sapply(inno_cor$inno1,FUN = function(x){unique(inno$midy[which(x == inno$q)])})
  inno_cor$inno2_midy = sapply(inno_cor$inno2,FUN = function(x){unique(inno$midy[which(x == inno$q)])})
  # set end and starting points for arrows
  inno_cor$x0 = unique(inno$midx) + 0.7*inno_radx
  inno_cor$x1 = unique(inno$midx) + 0.7*inno_radx
  inno_cor$y0 = ifelse(inno_cor$inno1 > inno_cor$inno2,inno_cor$inno1_midy + 0.75*inno_radx/asp, inno_cor$inno1_midy - 0.75*inno_radx/asp)
  inno_cor$y1 = ifelse(inno_cor$inno1 > inno_cor$inno2,inno_cor$inno2_midy - 0.75*inno_radx/asp, inno_cor$inno2_midy + 0.75*inno_radx/asp)

}

# for(i in 1:nrow(inno_cor)){
# diagram::curvedarrow(arr.type ="triangle", arr.length = arrHead_w,
#   curve = -1*inno_cor$curve[i], segment = c(0.1,0.9), arr.pos = 0.9, endhead = T, lwd = 0.7,
#   from = c(inno_cor$x0[i],inno_cor$y0[i]),
#   to = c(inno_cor$x1[i],inno_cor$y1[i]))
# diagram::curvedarrow(
#   arr.type ="triangle", arr.length = arrHead_w,
#   curve = inno_cor$curve[i], segment = c(0.1,0.9), arr.pos = 0.9, endhead = T, lwd = 0.7,
#   to = c(inno_cor$x0[i],inno_cor$y0[i]),
#   from = c(inno_cor$x1[i],inno_cor$y1[i]))
# }

# NODES FOR INTERACTION EFFECTS ================================================
if(infos$n_int > 0){
  int_pars <- infos$fix_pars_dyn[infos$fix_pars_dyn$isINT==1,]
  for(i in 1:nrow(int_pars)){
    # get positions of involved constructs
    preds <- w_nodes[(w_nodes$lag == int_pars$Lag[i] & w_nodes$construct == int_pars$Dpred[i]) |
                     (w_nodes$lag == int_pars$Lag2[i] & w_nodes$construct == int_pars$Dpred2[i]),]
    out.node <- w_nodes[w_nodes$lag == 0 & w_nodes$construct == int_pars$Dout[i],]

    # position offsets to avoid overlapping content
    int_x_offset = ifelse(int_pars$Lag[i] == int_pars$Lag2[i], 1.5*w_radx, 0)
    int_y_offset = ifelse(abs(diff(c(as.integer(int_pars$Dpred)[i], as.integer(int_pars$Lag2)[i]))) > 1, 1.5*w_radx/asp, 0)
    int_pars$midx[i] = get_mid_points(1, lims = preds$midx) + int_x_offset
    if(any(preds$lag == 0)){int_pars$midx[i] = preds$midx[preds$lag==0] - 1.5*w_radx}
    int_pars$midy[i] = get_mid_points(1, lims = preds$midy) + int_y_offset

    # paths to int node
    int_pars$p1_x0[i] <- preds$midx[1] + ifelse(preds$lag[1] > 0, w_radx, 0)
    int_pars$p1_y0[i] <- preds$midy[1]
    int_pars$p1_x1[i] <- int_pars$midx[i]
    int_pars$p1_y1[i] <- int_pars$midy[i]
    int_pars$p2_x0[i] <- preds$midx[2] + ifelse(preds$lag[2] > 0, w_radx, 0)
    int_pars$p2_y0[i] <- preds$midy[2]
    int_pars$p2_x1[i] <- int_pars$midx[i]
    int_pars$p2_y1[i] <- int_pars$midy[i]
    # path to dv node
    int_pars$x0[i] <- int_pars$midx[i]
    int_pars$y0[i] <- int_pars$midy[i]
    int_pars$x1[i] <- out.node$midx
    int_pars$y1[i] <- out.node$midy
    int_pars$arr_pos[i] <- get_arr_pos(int_pars$x0[i], int_pars$x1[i],int_pars$y0[i], int_pars$y1[i], radx = (w_radx+0.5*arrHead_w))
  }
}
# Arrows of lag 0 ==============================================================
w_arr <- infos$fix_pars_dyn[infos$fix_pars_dyn$Lag == "0",]
if(nrow(w_arr)>0){
  w_arr$Dpred = as.integer(w_arr$Dpred)
  w_arr$Dout = as.integer(w_arr$Dout)
  w_arr$Lag = as.integer(w_arr$Lag)
  for(i in 1:nrow(w_arr)){
    x_offset = ifelse(infos$maxLag > 1 | abs(diff(unlist(w_arr[i,c("Dout", "Dpred")])))>1, 0.2*w_radx,0)
    w_arr$x0[i] <- w_nodes[w_nodes$construct == w_arr$Dpred[i] & w_nodes$lag == (w_arr$Lag[i]), "midx"] + x_offset
    w_arr$x1[i] <- w_nodes[w_nodes$construct == w_arr$Dout[i] & w_nodes$lag == w_arr$Lag[i], "midx"] + x_offset
    w_arr$y0[i] <- w_nodes[w_nodes$construct == w_arr$Dpred[i] & w_nodes$lag == (w_arr$Lag[i]), "midy"]
    w_arr$y1[i] <- w_nodes[w_nodes$construct == w_arr$Dout[i] & w_nodes$lag == w_arr$Lag[i], "midy"]
    if(w_arr$Dpred[i] < w_arr$Dout[i]){
      w_arr$y0[i] <- w_arr$y0[i] - w_radx/asp/2
      w_arr$y1[i] <- w_arr$y1[i] + w_radx/asp
    }
    if(w_arr$Dpred[i] > w_arr$Dout[i]){
      w_arr$y0[i] <- w_arr$y0[i] + w_radx/asp/2
      w_arr$y1[i] <- w_arr$y1[i] - w_radx/asp
    }
  }
  w_arr_l0 <- w_arr
} else {
  w_arr_l0 <- NULL
}
# Arrows of lag 1 =============================================================

w_arr <- infos$fix_pars_dyn[infos$fix_pars_dyn$Lag == "1",]
w_arr$Dpred = as.integer(w_arr$Dpred)
w_arr$Dout = as.integer(w_arr$Dout)
w_arr$Lag = as.integer(w_arr$Lag)
for(i in 1:nrow(w_arr)){
  w_arr$x0[i] <- w_nodes[w_nodes$construct == w_arr$Dpred[i] & w_nodes$lag == (w_arr$Lag[i]), "midx"] + w_radx
  w_arr$x1[i] <- w_nodes[w_nodes$construct == w_arr$Dout[i] & w_nodes$lag+1 == w_arr$Lag[i], "midx"]
  w_arr$y0[i] <- w_nodes[w_nodes$construct == w_arr$Dpred[i] & w_nodes$lag == (w_arr$Lag[i]), "midy"]
  w_arr$y1[i] <- w_nodes[w_nodes$construct == w_arr$Dout[i] & w_nodes$lag+1 == w_arr$Lag[i], "midy"]
  w_arr$arr_pos[i] <- get_arr_pos(w_arr$x0[i], w_arr$x1[i],w_arr$y0[i], w_arr$y1[i], radx = (w_radx+0.5*arrHead_w))
}
w_arr_l1 <- w_arr
#### Greyed version if max_lag > 1
if(infos$maxLag > 1){
  w_arr <- infos$fix_pars_dyn[infos$fix_pars_dyn$Lag == "2",]
  w_arr$Dpred = as.integer(w_arr$Dpred)
  w_arr$Dout = as.integer(w_arr$Dout)
  w_arr$Lag = as.integer(w_arr$Lag)
  for(i in 1:nrow(w_arr)){
    w_arr$x0[i] <- w_nodes[w_nodes$construct == w_arr$Dpred[i] & w_nodes$lag == (w_arr$Lag[i]), "midx"] + w_radx
    w_arr$x1[i] <- w_nodes[w_nodes$construct == w_arr$Dout[i] & w_nodes$lag+1 == w_arr$Lag[i], "midx"]
    w_arr$y0[i] <- w_nodes[w_nodes$construct == w_arr$Dpred[i] & w_nodes$lag == (w_arr$Lag[i]), "midy"]
    w_arr$y1[i] <- w_nodes[w_nodes$construct == w_arr$Dout[i] & w_nodes$lag+1 == w_arr$Lag[i], "midy"]
    w_arr$arr_pos[i] <- get_arr_pos(w_arr$x0[i], w_arr$x1[i],w_arr$y0[i], w_arr$y1[i], radx = (w_radx+0.5*arrHead_w))
  }
  w_arr_l1.1 <- w_arr
}
##### HIGHER-ORDER LAGS ====
###### ARs
w_arr <- infos$fix_pars_dyn[infos$fix_pars_dyn$Lag == "2" & infos$fix_pars_dyn$isAR == 1,]
n_w_l2AR = nrow(w_arr)
if(nrow(w_arr) > 0){
  w_arr$Dpred = as.integer(w_arr$Dpred)
  w_arr$Dout = as.integer(w_arr$Dout)
  w_arr$Lag = as.integer(w_arr$Lag)
  for(i in 1:nrow(w_arr)){
    # horizontal lines
    w_arr$hline_x0[i] <- unique(w_nodes[w_nodes$construct == w_arr$Dpred[i] & w_nodes$lag == (w_arr$Lag[i]), "midx"])
    w_arr$hline_x1[i] <- unique(w_nodes[w_nodes$construct == w_arr$Dout[i] & w_nodes$lag+2 == (w_arr$Lag[i]), "midx"])
    y_adjust <- ifelse(w_arr$Dout[i] == max(w_arr$Dout),-w_radx/asp, w_radx/asp)
    w_arr$hline_y0[i] <- w_ys[w_arr$Dpred[i]] + 1.75*y_adjust
    w_arr$hline_y1[i] <- w_ys[w_arr$Dpred[i]] + 1.75*y_adjust
    # down arrows
    ## pred
    w_arr$pred_x0[i] <- w_arr$hline_x0[i]
    w_arr$pred_x1[i] <- w_arr$hline_x0[i]
    w_arr$pred_y0[i] <- w_ys[w_arr$Dpred[i]] + y_adjust
    w_arr$pred_y1[i] <- w_arr$hline_y1[i]
    ## pred
    w_arr$dv_x0[i] <- w_arr$hline_x1[i]
    w_arr$dv_x1[i] <- w_arr$hline_x1[i]
    w_arr$dv_y0[i] <- w_arr$hline_y0[i]
    w_arr$dv_y1[i] <- w_ys[w_arr$Dpred[i]] + y_adjust
  }
}
w_arr_l2.ar <- w_arr


# BETWEEN ======================================================================
n_nod_b <- infos$n_random
has_out <- ifelse(infos$n_out>0,1,0)
has_cov <- ifelse(infos$n_cov>1,1,0)
n_cov <- infos$n_cov-1
n_out <- infos$n_out
n_cov_paths <- nrow(infos$RE.PREDS)
n_out_paths <- nrow(infos$OUT)


# calculate positions
if(b_style == "h"){
  # positions of random effect par nodes
  b_poses_x = get_mid_points(n_nod_b, lims = c(begin.btw.x, end.btw.x))
  # basic settings (that may be overwritten based on other consitions)
  b_poses_y <- get_mid_points(8, lims = c(end.btw.y, begin.btw.y))[5]
  b_r_l_y   = get_mid_points(8, lims = c(end.btw.y, begin.btw.y))[3]
  b_cor_pos <- "t"
  # other vars ==============================================================
  cov_midy <- get_mid_points(4, lims = c(end.btw.y, begin.btw.y))[1]
  cov_midx <- get_mid_points(n_cov, lims = c(begin.btw.x, end.btw.x))
  out_midy <- get_mid_points(8, lims = c(end.btw.y, begin.btw.y))[8]
  out_midx <- get_mid_points(n_out, lims = c(begin.btw.x, end.btw.x))
  #### conditional
  if(has_out == 1 & has_cov == 0){
    b_poses_y <- get_mid_points(8, lims = c(end.btw.y, begin.btw.y))[4]
    b_r_l_y   = get_mid_points(8, lims = c(end.btw.y, begin.btw.y))[2]
    b_cor_pos <- "t"
  }
  if(has_out == 0 & has_cov == 1){
    cov_midy <- get_mid_points(10, lims = c(end.btw.y, begin.btw.y))[3]
    b_poses_y <- get_mid_points(10, lims = c(end.btw.y, begin.btw.y))[8]
    b_r_l_y   = get_mid_points(10, lims = c(end.btw.y, begin.btw.y))[10]
    b_cor_pos <- "b"
  }
  if(has_out == 1 & has_cov == 1){
    b_poses_y <- get_mid_points(9, lims = c(end.btw.y, begin.btw.y))[5]
    b_r_l_y   = get_mid_points(8, lims = c(end.btw.y, begin.btw.y))[6]
    b_cor_pos <- "b"
  }
  b_radx <- diff(b_poses_x[1:2])*scale_between
  ##### placement of correlations ==========================================

  b_cor_offset <- ifelse(has_cov == 1 & has_out == 1, b_radx, 0)
  b_cor_col <- ifelse(has_cov == 1 & has_out == 1, "grey", "black")
  b_cor_lty <- ifelse(has_cov == 1 & has_out == 1, "dashed", "solid")
  if(b_cor_pos == "t"){
    b_r_l_x <- c(min(b_poses_x) + b_cor_offset, max(b_poses_x) + b_cor_offset)
    b_r_arrows_x0 <- b_poses_x + b_cor_offset
    b_r_arrows_y0 <- rep(b_r_l_y,n_nod_b)
    b_r_arrows_x1 <- b_poses_x + 0.5*b_cor_offset
    b_r_arrows_y1 <- rep(b_poses_y,n_nod_b) + b_radx/asp
  }
  if(b_cor_pos == "b"){
    b_r_l_x <- c(min(b_poses_x)- b_cor_offset, max(b_poses_x)- b_cor_offset)
    b_r_arrows_x0 <- b_poses_x - b_cor_offset
    b_r_arrows_y0 <- rep(b_r_l_y,n_nod_b)
    b_r_arrows_x1 <- b_poses_x -b_cor_offset/2
    b_r_arrows_y1 <- rep(b_poses_y,n_nod_b) - b_radx/asp
  }
} else {        # vertical ....
  b_poses_x = get_mid_points(1, lims = c(end.btw.x, begin.btw.x))
  b_poses_y = get_mid_points(n_nod_b, lims = c(end.btw.y, begin.btw.y))
  b_radx <- diff(b_poses_y[1:2])*0.45
}
# get random pars -------
b_nodes = infos$re_pars
b_nodes$midx = b_poses_x
b_nodes$midy = b_poses_y
b_nodes$lab = plotmath_labeller(x = b_nodes$Param, y_fac_labs = y_fac_labs,
                                remove_lag_lab = remove_lag_lab, y_fac_lab_sep = y_fac_lab_sep)
b_nodes$radx = b_radx
if(has_cov == 1){
  covs <- infos$RE.PREDS[!(duplicated(infos$RE.PREDS$re_preds)),]
  covs$midx <- cov_midx
  covs$midy <- cov_midy
  covs$lab <- covs$re_preds
  arr.cov <- infos$RE.PREDS
  for(i in 1:nrow(arr.cov)){
    arr.cov$x0[i] <- unique(covs$midx[covs$re_preds == arr.cov$re_preds[i]])
    arr.cov$x1[i] <- unique(b_nodes$midx[b_nodes$Param == arr.cov$re_as_dv[i]])
    arr.cov$y0[i] <- unique(covs$midy[covs$re_preds == arr.cov$re_preds[i]]) - b_radx/asp
    arr.cov$y1[i] <- unique(b_nodes$midy[b_nodes$Param == arr.cov$re_as_dv[i]]) + b_radx/asp
  }
}
if(has_out == 1){
  outs <- infos$OUT[!(duplicated(infos$OUT$Var)),]
  outs$midx <- out_midx
  outs$midy <- out_midy
  outs$lab <- outs$Var
  arr.out <- infos$OUT
  for(i in 1:nrow(arr.out)){
    arr.out$x1[i] <- unique(outs$midx[outs$Var == arr.out$Var[i]])
    arr.out$x0[i] <- unique(b_nodes$midx[b_nodes$Param == arr.out$Pred[i]])
    arr.out$y1[i] <- unique(outs$midy[outs$Var == arr.out$Var[i]]) + b_radx/asp
    arr.out$y0[i] <- unique(b_nodes$midy[b_nodes$Param == arr.out$Pred[i]]) - b_radx/asp
  }
}
##### PLOTTING =================================================================
# add the midpoints
radx = 0.1 * diff(c(begin.decomp.x, end.decomp.x))
rady = radx * diff(fig_margins.x)/diff(fig_margins.y)

if(!is.null(file)){
  if(endsWith(file, suffix = ".png")){
    grDevices::png(filename = file, width = width, height = height, units = units, res = res, pointsize = pointsize, type = type, ...)
  }
}


diagram::openplotmat(graphics::par(mai = c(0,0,0,0)), ylim = fig_margins.y, xlim = fig_margins.x)
##### Separation lines :
diagram::straightarrow(from = c(end.decomp.x,fig_margins.y[1]), to = c(end.decomp.x, fig_margins.y[2]), lty="dashed", segment=c(0.05, 0.95),arr.type = "none", lwd=.9)
diagram::straightarrow(from = c(end.decomp.x,begin.wth.y),      to = c(end.wth.x, begin.wth.y),         lty="dashed", segment=c(   0, 0.95),arr.type = "none", lwd=.9)
##### Label Sections :
diagram::textempty(lab = "Between-Level", mid = c(end.decomp.x+w_radx*0.5,           end.btw.y-(radx/asp)), adj=0, font=2)
diagram::textempty(lab = "Within-Level",  mid = c(end.decomp.x+w_radx*0.5,           end.btw.y+(radx/asp)), adj=0, font=2)
diagram::textempty(lab = "Decomposition", mid = c(get_mid_points(1       , c(begin.decomp.x,end.decomp.x)), 0.9*fig_margins.y[2]), font=2)


##### DECOMPOSITION :

for ( i in 1:nrow(ind_nodes) ) {
  # manifest indicators
  diagram::textrect(shadow.size = 0,lwd=lwd_nodes, family=family,
    lab = parse(text=ind_nodes$lab[i]),
    mid = c(ind_nodes$midx[i], ind_nodes$midy[i]),
    radx = ind_nodes$ind_radx[i],
    rady = ind_nodes$ind_radx[i]/asp)

  if ( ind_nodes$is_wcen[i] == 1 ){

    # within-level latent factor(s)
    diagram::textellipse(shadow.size = 0,lwd=lwd_nodes, family=family,
      lab = parse(text=ind_nodes$Wf_lab[i]),
      mid = c(ind_nodes$Wf_midx[i], ind_nodes$Wf_midy[i]),
      radx = ind_nodes$indf_radx[i],
      rady = ind_nodes$indf_radx[i]/asp)

    # between-level latent factor(s)
    diagram::textellipse(shadow.size = 0,lwd=lwd_nodes, family=family,
      lab = parse(text=ind_nodes$Bf_lab[i]),
      mid = c(ind_nodes$Bf_midx[i], ind_nodes$Bf_midy[i]),
      radx = ind_nodes$indf_radx[i],
      rady = ind_nodes$indf_radx[i]/asp)

    # down arrows
    shape::Arrows(lwd = 0.9, lty = "solid", code = 2, arr.adj = 1, arr.type = "triangle", arr.length = arrHead_w,
      x0 = ind_nodes$d_arr_x0[i],
      y0 = ind_nodes$d_arr_y0[i],
      x1 = ind_nodes$d_arr_x1[i],
      y1 = ind_nodes$d_arr_y1[i])

    # up arrows
    shape::Arrows(lwd = 0.9, lty = "solid",code = 2, arr.adj = 1, arr.type = "triangle", arr.length = arrHead_w,
      x0 = ind_nodes$u_arr_x0[i],
      y0 = ind_nodes$u_arr_y0[i],
      x1 = ind_nodes$u_arr_x1[i],
      y1 = ind_nodes$u_arr_y1[i])

    # residuals and fixed loading labels
    if( infos$isLatent == 1 ){

      if( ind_nodes$sigmaW_isFree[i] == 1 ){
        diagram::straightarrow(lwd = 0.7, lty = "solid", arr.type ="triangle", endhead = T,
          from = c(ind_nodes$resW_arr_x0[i],ind_nodes$resW_arr_y0[i]),
          to   = c(ind_nodes$resW_arr_x1[i],ind_nodes$resW_arr_y1[i]),
          arr.pos = (1-0.8*arrHead_w), arr.length = 0.8*arrHead_w)
      }

      if( ind_nodes$sigmaB_isFree[i] == 1 ){
      diagram::straightarrow(lwd = 0.7, lty = "solid", arr.type ="triangle", endhead = T,
        from = c(ind_nodes$resB_arr_x0[i],ind_nodes$resB_arr_y0[i]),
        to   = c(ind_nodes$resB_arr_x1[i],ind_nodes$resB_arr_y1[i]),
        arr.pos = (1-0.8*arrHead_w), arr.length = 0.8*arrHead_w)
      }

      if( ind_nodes$lambdaW_isEqual[i] == "= 1" ){
        diagram::textempty(family=family,
          mid = c(ind_nodes$lamW1_midx[i],ind_nodes$lamW1_midy[i]),
          lab = ind_nodes$lamW1_lab[i], box.col = "white", cex = cex_loads, adj = adj_load_x
          )
      }

      if( ind_nodes$lambdaB_isEqual[i] == "= 1" ){
        diagram::textempty(family=family,
          mid = c(ind_nodes$lamB1_midx[i],ind_nodes$lamB1_midy[i]),
          lab = ind_nodes$lamB1_lab[i], box.col = "white", cex = cex_loads, adj = adj_load_x
        )
        }
      }
    }
  }

# END DECOMPOSITION -------------------------------------------------------------

#### WITHIN :


### lag 0 arrows
if(!is.null(w_arr_l0)){
  for(i in 1:nrow(w_arr_l0)){
    shape::Arrows(w_arr_l0$x0[i],w_arr_l0$y0[i], w_arr_l0$x1[i],w_arr_l0$y1[i],lwd = 0.9, lty = "solid", arr.type ="triangle", arr.adj = 1, arr.length = arrHead_w)
    # repeat arrows at each t
    x_diff = unique(w_nodes$midx[w_nodes$lag==0]) - unique(w_nodes$midx[w_nodes$lag==1])
    shape::Arrows(w_arr_l0$x0[i]-x_diff,w_arr_l0$y0[i], w_arr_l0$x0[i]-x_diff,w_arr_l0$y1[i],lwd = 0.9, lty = "solid", arr.type ="triangle", arr.adj = 1, arr.length = arrHead_w)
    if(w_arr_l0$isRandom[i]==1){
      diagram::straightarrow(from = c(w_arr_l0$x0[i],w_arr_l0$y0[i]), to = c(w_arr_l0$x1[i],w_arr_l0$y1[i]),lwd = 0.9, lty = "solid", arr.type ="circle", endhead = T, arr.pos = rand_dot_pos, arr.length = arrHead_w/2)
      diagram::straightarrow(from = c(w_arr_l0$x0[i]-x_diff,w_arr_l0$y0[i]), to = c(w_arr_l0$x1[i]-x_diff,w_arr_l0$y1[i]),lwd = 0.9, lty = "solid", arr.type ="circle", endhead = T, arr.pos = rand_dot_pos, arr.length = arrHead_w/2)
    }

  }
}
# interaction paths
if(infos$n_int>0){
  for(i in 1:infos$n_int){
    shape::Arrows(int_pars$p1_x0[i],int_pars$p1_y0[i], int_pars$p1_x1[i],int_pars$p1_y1[i],lwd = 0.9, lty = "solid", arr.adj = 1, arr.length =0)
    shape::Arrows(int_pars$p2_x0[i],int_pars$p2_y0[i], int_pars$p2_x1[i],int_pars$p2_y1[i],lwd = 0.9, lty = "solid", arr.adj = 1, arr.length =0)
    diagram::straightarrow(c(int_pars$x0[i],int_pars$y0[i]), c(int_pars$x1[i],int_pars$y1[i]),lwd = 0.9, lty = "solid", endhead = T, arr.length = arrHead_w, arr.pos = int_pars$arr_pos[i], arr.type = "triangle", arr.adj = 0.5)
    if(int_pars$isRandom[i]){
      diagram::straightarrow(from = c(int_pars$x0[i],int_pars$y0[i]), to = c(int_pars$x1[i],int_pars$y1[i]),lwd = 0.9, lty = "solid", arr.type ="circle", endhead = T, arr.pos = rand_dot_pos, arr.length = arrHead_w/2)
    }
  }
}

# innovation covariance
if(infos$n_inno_covs == 1){
  for(i in 1:nrow(inno_cov)){
    # add the factor loads
    diagram::textempty(lab = parse(text=inno_cov$loads[i]), mid = c(get_mid_points(1,c(inno_cov$midx[1], inno$x1[i])), get_mid_points(1,c(inno_cov$midy[1], inno$midy[i]))), lwd=lwd_nodes, family=family, adj = 0.5, cex = 0.8, box.col = "transparent")

    diagram::straightarrow(c(inno_cov$x0[i],inno_cov$y0[i]),c(inno_cov$x1[i],inno_cov$y1[i]), lwd=0.9,lty="solid",arr.pos = inno_cov$arr_pos[i],endhead=T,arr.length = arrHead_w,arr.type = "triangle", arr.adj = 1)
  }
  diagram::textellipse(lab = parse(text=inno_cov$lab[1]), mid = c(inno_cov$midx[1], inno_cov$midy[1]), radx = inno_radx, rady = inno_radx/asp, shadow.size = 0, lwd=lwd_nodes, family=family)
  # self-arrow
  diagram::curvedarrow(curve = 0.5,  from = c(inno_cov$x0[1]+3*inno_radx,inno_cov$y0[1]), to = c(inno_cov$x0[1]+1.5*inno_radx,inno_cov$y0[1]),lwd = 0.9, lty = "solid",endhead = T, arr.pos = 0.6, arr.length = arrHead_w, arr.type ="triangle")
  diagram::curvedarrow(curve = -0.5, from = c(inno_cov$x0[1]+3*inno_radx,inno_cov$y0[1]), to = c(inno_cov$x0[1]+1.5*inno_radx,inno_cov$y0[1]),lwd = 0.9, lty = "solid",endhead = T, arr.pos = 0.6, arr.length = arrHead_w, arr.type ="triangle")
  diagram::straightarrow(from = c(inno_cov$x0[1]+3*inno_radx,inno_cov$y0[1]), to = c(inno_cov$x1[1],inno_cov$y1[1]),lwd = 0.9, lty = "solid",  arr.type ="circle", endhead = T, arr.pos = 0, arr.length = arrHead_w/2)
}

# innovation correlations
if(infos$n_inno_cors > 0){
  for(i in 1:nrow(inno_cor)){
  diagram::curvedarrow(arr.type ="triangle", arr.length = arrHead_w,
                       curve = -1*inno_cor$curve[i], segment = c(0.1,0.9), arr.pos = 0.9, endhead = T, lwd = 0.7,
                       from = c(inno_cor$x0[i],inno_cor$y0[i]),
                       to = c(inno_cor$x1[i],inno_cor$y1[i]))
  diagram::curvedarrow(
    arr.type ="triangle", arr.length = arrHead_w,
    curve = inno_cor$curve[i], segment = c(0.1,0.9), arr.pos = 0.9, endhead = T, lwd = 0.7,
    to = c(inno_cor$x0[i],inno_cor$y0[i]),
    from = c(inno_cor$x1[i],inno_cor$y1[i]))
    }
  }


for(i in 1:nrow(w_nodes)){
  # within-level latent factor(s)
  if(w_nodes$w_cen[i] == 1){
    diagram::textellipse(lab = parse(text = w_nodes$lab[i]), mid = c(w_nodes$midx[i], w_nodes$midy[i]), radx = w_radx, rady = w_radx/asp, shadow.size = 0, lwd=lwd_nodes, family=family)
  } else {
    diagram::textrect(lab = parse(text = w_nodes$lab[i]), mid = c(w_nodes$midx[i], w_nodes$midy[i]), radx = w_radx*0.9, rady = w_radx*0.9/asp, shadow.size = 0, lwd=lwd_nodes, family=family)
  }
}
for(i in 1:nrow(inno)){
  # innovations
    diagram::textellipse(lab = parse(text=inno$lab[i]), mid = c(inno$midx[i], inno$midy[i]), radx = inno_radx, rady = inno_radx/asp, shadow.size = 0, lwd=lwd_nodes, family=family)
    # arrows
    shape::Arrows(x0=inno$x0[i],y0=inno$y0[i], x1=inno$x1[i],y1=inno$y1[i],lwd = 0.9, lty = "solid", arr.adj = 1, arr.length = arrHead_w, arr.type = "triangle")
    # self-arrow
    diagram::curvedarrow(curve = 0.5,  from = c(inno$x0[i]+3*inno_radx,inno$y0[i]), to = c(inno$x0[i]+1.5*inno_radx,inno$y0[i]),lwd = 0.9, lty = "solid",endhead = T, arr.pos = 0.6, arr.length = arrHead_w, arr.type ="triangle")
    diagram::curvedarrow(curve = -0.5, from = c(inno$x0[i]+3*inno_radx,inno$y0[i]), to = c(inno$x0[i]+1.5*inno_radx,inno$y0[i]),lwd = 0.9, lty = "solid",endhead = T, arr.pos = 0.6, arr.length = arrHead_w, arr.type ="triangle")
    if(inno$isRandom[i]==1){
      diagram::straightarrow(from = c(inno$x0[i]+3*inno_radx,inno$y0[i]), to = c(inno$x1[i],inno$y1[i]),lwd = 0.9, lty = "solid",  arr.type ="circle", endhead = T, arr.pos = 0, arr.length = arrHead_w/2)
    }
}
### lag 1 arrows
for(i in 1:nrow(w_arr_l1)){
  if(w_arr_l1$isAR[i] == 1 | w_arr_l1$y0[i] == w_arr_l1$y1[i]){
    shape::Arrows(w_arr_l1$x0[i],w_arr_l1$y0[i], (w_arr_l1$x1[i]-w_radx), w_arr_l1$y1[i], lwd = 0.9, lty = "solid", arr.type ="triangle", arr.length = arrHead_w, arr.adj = 1)
  } else {
    diagram::straightarrow(c(w_arr_l1$x0[i],w_arr_l1$y0[i]), c(w_arr_l1$x1[i], w_arr_l1$y1[i]), lwd = 0.9, lty = "solid", endhead = T, arr.type ="triangle", arr.pos = w_arr_l1$arr_pos[i], arr.length = arrHead_w, arr.adj = 0.5)
  }
  if(w_arr_l1$isRandom[i]==1){
    diagram::straightarrow(from = c(w_arr_l1$x0[i],w_arr_l1$y0[i]), to = c(w_arr_l1$x1[i],w_arr_l1$y1[i]),lwd = 0.9, lty = "solid", arr.type ="circle", endhead = T, arr.pos = rand_dot_pos, arr.length = arrHead_w/2)
  }
  if(infos$maxLag > 1){
    if(w_arr_l1$isAR[i] == 1){
      shape::Arrows(w_arr_l1.1$x0[i],w_arr_l1.1$y0[i], (w_arr_l1.1$x1[i]-w_radx), w_arr_l1.1$y1[i], lwd = 0.9, lcol = "grey", lty = "solid", arr.type ="triangle", arr.length = arrHead_w, arr.adj = 1)
    } else {
      diagram::straightarrow(c(w_arr_l1.1$x0[i],w_arr_l1.1$y0[i]), c(w_arr_l1.1$x1[i], w_arr_l1.1$y1[i]), lwd = 0.9, lty = "solid", lcol = "grey", endhead = T, arr.type ="triangle", arr.pos = w_arr_l1.1$arr_pos[i], arr.length = arrHead_w, arr.adj = 0.5)
    }
    if(w_arr_l1$isRandom[i]==1){
      diagram::straightarrow(from = c(w_arr_l1.1$x0[i],w_arr_l1.1$y0[i]), to = c(w_arr_l1.1$x1[i],w_arr_l1.1$y1[i]),lwd = 0.8, lty = "solid", lcol = "grey", arr.type ="circle", arr.col = "grey", endhead = T, arr.pos = rand_dot_pos, arr.length = arrHead_w/2)
    }
  }
}
if(n_w_l2AR > 0){
  for(i in 1:nrow(w_arr_l2.ar)){
    diagram::straightarrow(from = c(w_arr_l2.ar$hline_x0[i],w_arr_l2.ar$hline_y0[i]), to = c(w_arr_l2.ar$hline_x1[i],w_arr_l2.ar$hline_y1[i]),lwd = 0.9, lty = "solid", arr.type ="none")
    if(w_arr_l2.ar$isRandom[i] == 1){
      diagram::straightarrow(from = c(w_arr_l2.ar$hline_x0[i],w_arr_l2.ar$hline_y0[i]), to = c(w_arr_l2.ar$hline_x1[i],w_arr_l2.ar$hline_y1[i]),lwd = 0.9, lty = "solid", arr.type ="circle",arr.pos = (1)/2, arr.length = arrHead_w/2)
    }
    diagram::straightarrow(from = c(w_arr_l2.ar$pred_x0[i],w_arr_l2.ar$pred_y0[i]), to = c(w_arr_l2.ar$pred_x1[i],w_arr_l2.ar$pred_y1[i]),lwd = 0.9, lty = "solid", arr.type ="none")
    diagram::straightarrow(from = c(w_arr_l2.ar$dv_x0[i],w_arr_l2.ar$dv_y0[i]), to = c(w_arr_l2.ar$dv_x1[i],w_arr_l2.ar$dv_y1[i]),lwd = 0.9, lty = "solid", arr.type ="triangle",endhead = T, arr.pos = (1-arrHead_w), arr.length = arrHead_w)
  }
}
# interactions
if(infos$n_int>0){
  for(i in 1:infos$n_int){
    diagram::textellipse(mid = c(int_pars$midx[i], int_pars$midy[i]), radx = scale_int*w_radx, rady =  scale_int*w_radx/asp, shadow.size = 0, lwd=lwd_nodes, family=family)
  }
}

#### BETWEEN :

for ( i in 1:n_nod_b ) {

  # manifest indicators
  diagram::textellipse(shadow.size = 0, cex=cex_b, lwd=lwd_nodes,family=family,
    lab = parse(text = b_nodes$lab[i]),
    mid = c(b_nodes$midx[i], b_nodes$midy[i]),
    radx = b_nodes$radx[i],
    rady = b_nodes$radx[i]/asp)

  # correlation "arrows"
  shape::Arrows(
    b_r_l_x[1],b_r_l_y[1], b_r_l_x[2],b_r_l_y[1],
    lwd = 0.9, lty = b_cor_lty, lcol = b_cor_col, arr.length = 0)
  shape::Arrows(
    b_r_arrows_x0[i],b_r_arrows_y0[i], b_r_arrows_x1[i],b_r_arrows_y1[i],
    lwd = 0.9, lty = b_cor_lty, arr.type ="triangle",
    arr.adj = 1, arr.length = arrHead_b, lcol = b_cor_col)
}
# re predictors
if(has_cov == 1){
  for(i in 1:n_cov){
    diagram::textrect(lab = covs$lab[i],mid = c(covs$midx[i], covs$midy[i]), radx = b_nodes$radx[i], rady = b_nodes$radx[i]/asp, shadow.size = 0, cex=cex_b, lwd=lwd_nodes,family=family)
  }
  for(i in 1:nrow(arr.cov)){
    diagram::straightarrow(from = c(arr.cov$x0[i],arr.cov$y0[i]), to = c(arr.cov$x1[i],arr.cov$y1[i]),lwd = 0.9, lty = "solid", arr.type ="triangle", endhead = T, arr.pos = 1-0.3*arrHead_b, arr.length = arrHead_b)
  }
}
# outcomes
if(has_out == 1){
  for(i in 1:n_out){
    diagram::textrect(lab = outs$lab[i],mid = c(outs$midx[i], outs$midy[i]), radx = b_nodes$radx[i], rady = b_nodes$radx[i]/asp, shadow.size = 0, cex=cex_b, lwd=lwd_nodes, family=family)
  }
  for(i in 1:nrow(arr.out)){
    diagram::straightarrow(from = c(arr.out$x0[i],arr.out$y0[i]), to = c(arr.out$x1[i],arr.out$y1[i]),lwd = 0.9, lty = "solid", arr.type ="triangle", endhead = T, arr.pos = 1-0.75*arrHead_b, arr.length = arrHead_b)
  }
}

if(!is.null(file)){
  grDevices::dev.off()
}

}

