#' Create TeX Model Formula from VARmodel object
#'
#' @param VARmodel The VARmodel object.
#' @param data An optional `data.frame` including the variables to be used
#' in the formula.
#' @param labels An optional character string including the names to be used
#' in the formula.
#'
#' @return An Rmarkdown file that is automatically rendered to a pdf document.
#' @export
#'
#' @examples
VARmodelformula <- function(VARmodel, data = NULL, labels = NULL) {

  # extract model infos
  infos <- VARmodelEval(VARmodel)

  # extract model data frame
  model <- VARmodel

  # create empty markdown file, delete if already existing
  if (file.exists("formula.rmd")) {
    file.remove("formula.rmd")
  }
  rmarkdown::draft(file = "formula.rmd",
                   template = "formula",
                   package = "mlts",
                   edit = FALSE)

  # latex center begin and end
  begin_center <- "\\begin{center}"
  end_center <- "\\end{center}"
  # latex math begin and end
  begin_math <- "\\[\n\\begin{aligned}"
  end_math <- "\\end{aligned}\n\\]"
  # latex bmatrix begin and end
  begin_bmatrix <- "\\begin{bmatrix}"
  end_bmatrix <- "\\end{bmatrix}"


  # decomposition #############################################################

  # caption
  dcf_caption <- "Decomposition."

  # for manifest time-series constructs
  if (infos$isLatent == F) {
    ts_vec <- c()
    between_vec <- c()
    within_vec <- c()
    for (i in 1:infos$q) {
      ts_vec <- c(ts_vec, paste0("y_{", i, ", t} \\\\"))
      between_vec <- c(between_vec, paste0("\\mu_{", i, "} \\\\"))
      within_vec <- c(within_vec, paste0("y_{", i, ", t}^w \\\\"))
    }

    # paste together
    ts <- paste(ts_vec, collapse = "\n")
    between <- paste(between_vec, collapse = "\n")
    within <- paste(within_vec, collapse = "\n")
    # formulae
    dcf_lhs <- paste(begin_bmatrix, ts, end_bmatrix, sep = "\n")
    dcf_rhs <- paste(
      begin_bmatrix, between, end_bmatrix, "+",
      begin_bmatrix, within, end_bmatrix,
      sep = "\n"
    )

    # decomposition formula
    dcf <- paste(begin_math, dcf_lhs, "&=", dcf_rhs, end_math)

  } else {
    # for latent time-series constructs
    ts_vec <- c()
    ts_w_vec <- c() # within-model indicators
    ts_b_vec <- c() # between-model indicators
    lat_w_vec <- c() # within-model latent variables
    lat_b_vec <- c() # between-model latent variables
    eps_w_vec <- c() # within-model residuals
    eps_b_vec <- c() # between-model residuals
    lam_w_list <- list()
    lam_b_list <- list()
    lam_w_mat <- matrix(0, ncol = infos$q, nrow = sum(infos$p)) # within-model loadings
    lam_b_mat <- matrix(ncol = infos$q, nrow = sum(infos$p)) # between-model loadings
    int_b_vec <- c() # between-model intercepts
    for (i in 1:infos$q) {
      lat_w_vec <- c(lat_w_vec, paste0("\\eta^w_{", i, ", t} \\\\"))
      lat_b_vec <- c(lat_b_vec, paste0("\\eta^b_{", i, "} \\\\"))
      lam_w_list[[i]] <- matrix(ncol = infos$q, nrow = infos$p[i], data = 0)
      lam_b_list[[i]] <- matrix(ncol = infos$q, nrow = infos$p[i], data = 0)
      for (j in 1:infos$p[i]) {
        # decomposition to within- and between-parts
        ts_vec <- c(ts_vec, paste0("y_{", i, j, ", t} \\\\"))
        ts_w_vec <- c(ts_w_vec, paste0("y_{", i, j, ", t}^w \\\\"))
        ts_b_vec <- c(ts_b_vec, paste0("\\mu_{", i, j, "} \\\\"))
        eps_w_vec <- c(eps_w_vec, paste0("\\varepsilon_{", i, j, ", t}^w \\\\"))
        eps_b_vec <- c(eps_b_vec, paste0("\\varepsilon_{", i, j, "}^b \\\\"))
        int_b_vec <- c(int_b_vec, paste0("\\alpha_{", i, j, "}^b \\\\"))
        # within-model loading matrix
        if (
          infos$indicators[
            infos$indicators$q == i & infos$indicators$p == j,
            "lambdaW_isFree"
          ] == 0
        ) {
          lam_w_list[[i]][j, i] <- "1"
        } else {
          lam_w_list[[i]][j, i] <- paste0("\\lambda_{", i, j, "}^w")
        }
        # between-model loading matrix
        if (
          infos$indicators[
            infos$indicators$q == i & infos$indicators$p == j,
            "lambdaB_isFree"
          ] == 0
        ) {
          lam_b_list[[i]][j, i] <- "1"
        } else {
          lam_b_list[[i]][j, i] <- paste0("\\lambda_{", i, j, "}^b")
        }
      }
    }
    # bind matrices
    lam_w_mat <- do.call(rbind, lam_w_list)
    lam_w_vec <- c()
    for (i in 1:nrow(lam_w_mat)) {
      for (j in 1:ncol(lam_w_mat)) {
        if (j == ncol(lam_w_mat)) {
          lam_w_vec <- c(lam_w_vec, paste0(lam_w_mat[i, j], " \\\\"))
        } else {
          lam_w_vec <- c(lam_w_vec, paste0(lam_w_mat[i, j], " &"))
        }
      }
    }
    lam_b_mat <- do.call(rbind, lam_b_list)
    lam_b_vec <- c()
    for (i in 1:nrow(lam_b_mat)) {
      for (j in 1:ncol(lam_b_mat)) {
        if (j == ncol(lam_b_mat)) {
          lam_b_vec <- c(lam_b_vec, paste0(lam_b_mat[i, j], " \\\\"))
        } else {
          lam_b_vec <- c(lam_b_vec, paste0(lam_b_mat[i, j], " &"))
        }
      }
    }

    # paste together
    ts_w <- paste(ts_w_vec, collapse = "\n")
    ts_b <- paste(ts_b_vec, collapse = "\n")
    lam_w <- paste(lam_w_vec, collapse = "\n")
    lam_b <- paste(lam_b_vec, collapse = "\n")
    int_b <- paste(int_b_vec, collapse = "\n")
    lat_w <- paste(lat_w_vec, collapse = "\n")
    lat_b <- paste(lat_b_vec, collapse = "\n")
    eps_w <- paste(eps_w_vec, collapse = "\n")
    eps_b <- paste(eps_b_vec, collapse = "\n")

    # formulae
    dcf_lhs <- paste(begin_bmatrix, ts, end_bmatrix, sep = "\n")
    dcf_rhs <- paste(
      begin_bmatrix, ts_b, end_bmatrix, "+",
      begin_bmatrix, ts_w, end_bmatrix,
      sep = "\n"
    )
    dcf_lhs_w <- paste(begin_bmatrix, ts_w, end_bmatrix, sep = "\n")
    dcf_lhs_b <- paste(begin_bmatrix, ts_b, end_bmatrix, sep = "\n")
    dcf_rhs_w <- paste(
      begin_bmatrix, lam_w, end_bmatrix,
      begin_bmatrix, lat_w, end_bmatrix, "+",
      begin_bmatrix, eps_w, end_bmatrix,
      sep = "\n"
    )
    dcf_rhs_b <- paste(
      begin_bmatrix, int_b, end_bmatrix, "+",
      begin_bmatrix, lam_b, end_bmatrix,
      begin_bmatrix, lat_b, end_bmatrix, "+",
      begin_bmatrix, eps_b, end_bmatrix,
      sep = "\n"
    )

    # decomposition formula with within- and between-level measurement model
    dcf <- paste(
      begin_math,
      dcf_lhs, "&=", dcf_rhs, "\\\\",
      dcf_lhs_w, "&=", dcf_rhs_w, "\\\\",
      dcf_lhs_b, "&=", dcf_rhs_b, "\\\\",
      end_math
    )
  }


  # with caption
  decomposition <- paste(begin_center, dcf_caption, end_center, dcf, sep = "\n")


  # paste within_model to markdown
  cat(decomposition, file = "formula.rmd", append = TRUE)



  # within-model ##############################################################

  # caption
  wmf_caption <- "Within-model."

  # initiate empty dv vector
  dvs_vec <- c()
  # left hand side
  for (i in 1:infos$q) {
    dvs_vec <- c(dvs_vec, paste0("y_{", i, ", t}^w \\\\"))
  }
  dvs <- paste(dvs_vec, collapse = "\n")
  wmf_lhs <- paste(begin_bmatrix, dvs, end_bmatrix, sep = "\n")

  # autoregressive / cross-lagged parameter matrix

  # store number of phi-parameters for loops
  all_phis <- model[grepl("phi", model$Param) & grepl("Fix", model$Type), ]

  phi_mat_vec <- c()
  # loop across phi subscripts
  for (i in 1:infos$q) {
    for (j in 1:infos$q) {
      # if phi_ij exists in model frame (i.e., is not fixed to zero), paste in matrix
      if (any(all_phis$Param == paste0("phi_", i, j))) {
        if (j == infos$q) {
          phi_mat_vec <- c(phi_mat_vec, paste0("\\phi_{", i, j, "} \\\\ \n"))
        } else if (i == VARmodel$q & j == VARmodel$q) {
          # line before end of bmatrix must not be broken
          phi_mat_vec <- c(phi_mat_vec, paste0("\\phi_{", i, j, "} \\\\"))
        } else {
          phi_mat_vec <- c(phi_mat_vec, paste0("\\phi_{", i, j, "} & "))
        }
      } else {
        # if phi_ij does not exist in model frame, paste empty cell
        if (j == infos$q) {
          phi_mat_vec <- c(phi_mat_vec, paste0("0 \n"))
        } else if (i == infos$q & j == infos$q) {
          # line before end of bmatrix must not be broken
          phi_mat_vec <- c(phi_mat_vec, paste0("0 \\\\"))
        } else {
          phi_mat_vec <- c(phi_mat_vec, paste0("0 &"))
        }
      }
    }
  }
  phi_mat <- paste(phi_mat_vec, collapse = "")

  # initiate empty time series vector
  ts_vec <- c()
  # time series loop
  for (i in 1:infos$q) {
    ts_vec <- c(ts_vec, paste0("y_{", i, ", t - 1}^w \\\\")
    )
  }
  ts <- paste(ts_vec, collapse = "\n")

  # initiate empty innovation vector
  innos_vec <- c()
  # innovations loop
  for (i in 1:infos$q) {
    innos_vec <- c(innos_vec, paste0("\\zeta_{y,", i, ", t} \\\\"))
  }
  innos <- paste(innos_vec, collapse = "\n")

  # create right hand side formula
  wmf_rhs <- paste(
    begin_bmatrix, phi_mat, end_bmatrix,
    begin_bmatrix, ts, end_bmatrix, "+",
    begin_bmatrix, innos, end_bmatrix,
    collapse = "\n"
  )

  inno_dist <- ",\\\\ \n & \\text{with}~
  \\zeta_{y, i} \\sim \\mathit{MVN}(\\mathbf{0}, \\mathbf{\\Psi})"

  # within-model formula
  wmf <- paste(begin_math, wmf_lhs, "&=", wmf_rhs, inno_dist, end_math)

  # with caption
  within_model <- paste(begin_center, wmf_caption, end_center, wmf)


  # paste within_model to markdown
  cat(within_model, file = "formula.rmd", append = TRUE)


  # between-model #############################################################

  # caption
  bmf_caption <- "Between-model."

  all_bpars <- model[grepl("Fix", model$Type), ]
  # delete this if parameter names with underscore are implemented
  # replaces dots with nothing
  all_bpars$Param <- gsub("\\.", "", all_bpars$Param)
  names_bpars <- gsub(
    # replace underscore digit with latex subscript
    "()_(\\d+)", "\\1_{\\2}", gsub(
      # replace lnsigma with log(pi)
      "(lnsigma2|lnsigma)_(\\d+)", "log(\\\\pi_{\\2})", all_bpars$Param
    )
  )
  # one version unit-specific
  names_bpars_i <- gsub(
    # replace underscore digit with latex subscript
    "()_(\\d+)", "\\1_{\\2,i}", gsub(
      # replace lnsigma with log(pi)
      "(lnsigma2|lnsigma)_(\\d+)", "log(\\\\pi_{\\2,i})", all_bpars$Param
    )
  )

  n_bpars <- nrow(all_bpars)
  # left-hand side
  bpars_vec <- c()
  for (i in 1:n_bpars) {
    bpars_vec <- c(bpars_vec, paste0("\\", names_bpars_i[i], "\\\\"))
  }
  bpars <- paste(bpars_vec, collapse = "\n")
  # create between model left hand side
  bmf_lhs <- paste(begin_bmatrix, bpars, end_bmatrix, collapse = "\n")

  # right-hand side fixed effects
  fixef_vec <- c()
  for (i in 1:n_bpars) {
    fixef_vec <- c(fixef_vec, paste0("\\gamma_{\\", names_bpars[i], "}\\\\"))
  }
  fixef <- paste(fixef_vec, collapse = "\n")

  # right-hand side random effects
  ranef_vec <- c()
  for (i in 1:n_bpars) {
    if (all_bpars$isRandom[i] == 1) {
      ranef_vec <- c(ranef_vec, paste0("\\upsilon_{\\", names_bpars[i], ",i}\\\\"))
    } else {
      ranef_vec <- c(ranef_vec, paste0("0\\\\"))
    }
  }
  ranef <- paste(ranef_vec, collapse = "\n")

  # create right hand side formula
  bmf_rhs <- paste(
    begin_bmatrix, fixef, end_bmatrix, "+",
    begin_bmatrix, ranef, end_bmatrix,
    collapse = "\n"
  )

  ranef_dist <- ",\\\\ \n & \\text{with}~
  \\upsilon_{i} \\sim \\mathit{MVN}(\\mathbf{0}, \\mathbf{\\Omega})"

  # between-model formula
  bmf <- paste(begin_math, bmf_lhs, "&=", bmf_rhs, ranef_dist, end_math)

  # with caption
  between_model <- paste(begin_center, bmf_caption, end_center, bmf)


  # paste within_model to markdown
  cat(between_model, file = "formula.rmd", append = TRUE)

  # render markdown input #####################################################

  # render markdown
  rmarkdown::render(
    input = "formula.rmd"
  )

}
