#' Create TeX Model Formula from mlts model object
#'
#' @param model A model built with \code{\link[mlts]{mlts_model}}.
#' @param file An optional string containing the name of the file and file path.
#' Has to end with .pdf file format.
#' @param keep_tex Logical. Should the TeX file be kept (additional to the
#' Rmd file)? Defaults to `FALSE`.
#' @param ts To be included in future releases.
#' An optional character vector containing the names of the time-series
#' variables or indicators.
#' @param covariates To be included in future releases.
#' An optional character vector containing the names of the between-level
#' covariates.
#' @param outcomes To be included in future releases.
#' An optional character vector containing the names of the between-level
#' outcomes.
#' @return An RMarkdown file that is automatically rendered to a pdf document.
#' @export
#'
#' @examples
#' \donttest{
#' # build a simple vector-autoregressive mlts model with two time-series variables
#' var_model <- mlts_model(q = 2)
#'
#' # create formula from the specified model
#' mlts_model_formula(model = var_model)
#' }
#'
mlts_model_formula <- function(model, file = NULL,
                               keep_tex = FALSE,
                               ts = NULL, covariates = NULL,
                               outcomes = NULL) {

  # extract model infos
  infos <- mlts_model_eval(model)

  # extract model data frame
  model <- model

  # get file path, choose top-level directory if not specified
  if (!is.null(file)) {
    rmd_file <- gsub(".pdf", ".rmd", file)
    pdf_file <- file
  } else {
    rmd_file <- "formula.rmd"
    pdf_file <- "formula.pdf"
  }

  # create empty markdown file, delete if already existing
  if (file.exists(rmd_file)) {
    file.remove(rmd_file)
  }
  rmarkdown::draft(file = rmd_file,
                   template = "formula",
                   package = "mlts",
                   edit = FALSE)

  # add keep_tex to YAML frontmatter if desired
  if (keep_tex == TRUE) {
    yaml_keep_tex <- gsub(
      "output: pdf_document", "output: pdf_document:\n\tkeep_tex: true",
      readLines(rmd_file)
    )
    writeLines(yaml_keep_tex, rmd_file)
  }

  # latex center begin and end
  begin_center <- "\\begin{center}"
  end_center <- "\\end{center}"
  # latex math begin and end
  begin_math <- "\\[\n\\begin{gathered}"
  end_math <- "\\end{gathered}\n\\]\n\\vspace{1em}\n"
  # latex bmatrix begin and end
  begin_bmatrix <- "\\begin{bmatrix}"
  end_bmatrix <- "\\end{bmatrix}"


  # decomposition #############################################################

  # caption
  dcf_caption <- "Decomposition."

  # for manifest time-series constructs
  if (infos$isLatent == FALSE) {
    ts_vec <- c()
    between_vec <- c()
    within_vec <- c()
    for (i in 1:infos$q) {
      ts_vec <- c(ts_vec, paste0("y_{", i, ", it} \\\\"))
      between_vec <- c(between_vec, paste0("\\mu_{", i, ",i} \\\\"))
      within_vec <- c(within_vec, paste0("y_{", i, ", it}^w \\\\"))
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
    dcf <- paste(begin_math, dcf_lhs, "=", dcf_rhs, end_math)

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
    lam_b_mat <- matrix(0, ncol = infos$q, nrow = sum(infos$p)) # between-model loadings
    int_b_vec <- c() # between-model intercepts
    for (i in 1:infos$q) {
      lat_w_vec <- c(lat_w_vec, paste0("\\eta^w_{", i, ", it} \\\\"))
      lam_w_list[[i]] <- matrix(ncol = infos$q, nrow = infos$p[i], data = 0)
      if (all(infos$indicators[infos$indicators$q == i, "btw_factor"] == 1)) {
        # if common between-factor is modeled, use eta
        lat_b_vec <- c(lat_b_vec, paste0("\\eta^b_{", i, ",i} \\\\"))
        # and use loading matrix
        lam_b_list[[i]] <- matrix(ncol = infos$q, nrow = infos$p[i], data = 0)
      } else {
        # if no common between-factor is modeled, use mu (nothing)
        lat_b_vec <- lat_b_vec
        # and don't use loading matrix
        lam_b_list[[i]] <- NULL
      }
      for (j in 1:infos$p[i]) {
        # decomposition to within- and between-parts
        ts_vec <- c(ts_vec, paste0("y_{", i, j, ", it} \\\\"))
        ts_w_vec <- c(ts_w_vec, paste0("y_{", i, j, ", it}^w \\\\"))
        ts_b_vec <- c(ts_b_vec, paste0("\\mu_{", i, j, ",i} \\\\"))
        # if infos$p[i] is 1, add 0 as within-level epsilon
        if (infos$p[i] == 1) {
          eps_w_vec <- c(eps_w_vec, paste0("0 \\\\"))
        } else {
          eps_w_vec <- c(eps_w_vec, paste0("\\varepsilon_{", i, j, ", it}^w \\\\"))
        }
        # if common between-factor is modeled, build between-level formula
        if (all(infos$indicators[infos$indicators$q == i, "btw_factor"] == 1)) {
          # if common between-factor is modeled, use eta
          int_b_vec <- c(int_b_vec, ifelse(
            j == 1,
            "0 \\\\", # first indicator intercept fixed to zero
            paste0("\\alpha_{", i, j, ",i}^b \\\\")
          ))
          # and residual vector
          eps_b_vec <- c(eps_b_vec, ifelse(
            j == 1,
            "0 \\\\", # first residual fixed to zero
            paste0("\\varepsilon_{", i, j, ",i}^b \\\\")
          ))
        } else {
          # if no common between-factor is modeled, don't add formula
          int_b_vec <- int_b_vec
          eps_b_vec <- eps_b_vec
        }
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
        if (all(infos$indicators[infos$indicators$q == i, "btw_factor"] == 1)) {
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
        } else {
          lam_b_list <- lam_b_list
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
    if (!is.null(lam_b_mat)) {
      for (i in 1:nrow(lam_b_mat)) {
        for (j in 1:ncol(lam_b_mat)) {
          if (j == ncol(lam_b_mat)) {
            lam_b_vec <- c(lam_b_vec, paste0(lam_b_mat[i, j], " \\\\"))
          } else {
            lam_b_vec <- c(lam_b_vec, paste0(lam_b_mat[i, j], " &"))
          }
        }
      }
    }

    # paste together
    ts <- paste(ts_vec, collapse = "\n")
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
    dcf_rhs_w <- paste(
      begin_bmatrix, lam_w, end_bmatrix,
      begin_bmatrix, lat_w, end_bmatrix, "+",
      begin_bmatrix, eps_w, end_bmatrix,
      sep = "\n"
    )
    if (!is.null(lat_b_vec)) {
      dcf_lhs_b <- paste(begin_bmatrix, ts_b, end_bmatrix, sep = "\n")
      dcf_rhs_b <- paste(
        begin_bmatrix, int_b, end_bmatrix, "+",
        begin_bmatrix, lam_b, end_bmatrix,
        begin_bmatrix, lat_b, end_bmatrix, "+",
        begin_bmatrix, eps_b, end_bmatrix,
        sep = "\n"
      )
    }

    # residual distribution
    if (any(infos$p > 1)) {
      eps_dist_w <- paste0(
        ",\\\\ \n\\text{with}~",
        "\\varepsilon^w_{jk, it} \\sim \\mathit{N}(0,",
        "\\sigma^2_{\\varepsilon^w_{jk}})~",
        "\\text{for indicator $j$ of construct $k$}\\\\"
      )
    } else {
      eps_dist_w <- "\\\\"
    }

    if (any(infos$p > 1)) {
      eps_dist_b <- paste0(
        ",\\\\ \n\\text{with}~",
        "\\varepsilon^b_{jk,i} \\sim \\mathit{N}(0,",
        "\\sigma^2_{\\varepsilon^b_{jk}})~",
        "\\text{for indicator $j$ of construct $k$}\\\\"
      )
    } else {
      eps_dist_b <- NULL
    }


    # decomposition formula with within- and between-level measurement model
    dcf <- paste(
      begin_math,
      dcf_lhs, "=", dcf_rhs, "\\\\",
      dcf_lhs_w, "=", dcf_rhs_w,
      # if (!is.null(eps_dist_w)) {
      #   paste(eps_dist_w)
      # },
      eps_dist_w,
      if (!is.null(lat_b_vec)) {
        if (!is.null(eps_dist_b)) {
          paste(dcf_lhs_b, "=", dcf_rhs_b, eps_dist_b)
        } else {
          paste(dcf_lhs_b, "=", dcf_rhs_b)
        }
      },
      end_math
    )
  }


  # with caption
  decomposition <- paste(begin_center, dcf_caption, end_center, dcf, sep = "\n")


  # paste within_model to markdown
  cat(decomposition, file = rmd_file, append = TRUE)



  # within-model ##############################################################

  # caption
  wmf_caption <- "Within-level model."

  # initiate empty dv vector
  dvs_vec <- c()
  # left hand side
  for (i in 1:infos$q) {
    dvs_vec <- c(dvs_vec, ifelse(
      infos$isLatent == TRUE,
      paste0("\\eta_{", i, ", it}^w \\\\"),
      paste0("y_{", i, ", it}^w \\\\")
    ))
  }
  dvs <- paste(dvs_vec, collapse = "\n")
  wmf_lhs <- paste(begin_bmatrix, dvs, end_bmatrix, sep = "\n")

  # autoregressive / cross-lagged parameter matrix

  # store number of phi-parameters for loops
  all_phis <- model[grepl("phi", model$Param) & grepl("Fix", model$Type), ]
  all_phis$names <- gsub(
    # replace underscore digit with latex subscript
    "(\\w+)(\\(\\d\\))_(\\d+)", "\\\\\\1_{\\2\\3}", all_phis$Param
  )

  # all_phis$mat_row <- as.numeric(
  #   gsub(
  #     "(\\w+)\\((\\d)\\)_(\\d)(\\d)", "\\3", all_phis$Param
  #   )
  # )
  # think about this if you have time to spare
  # how can matrix column be calculated from capture group
  # \\2 and \\3 alone?
  # all_phis$mat_col <- as.numeric(
  #   gsub(
  #     "(\\w+)\\((\\d)\\)_(\\d)(\\d)", "\\4", all_phis$Param
  #   )
  # ) ^ as.numeric(
  #   gsub(
  #     "(\\w+)\\((\\d)\\)_(\\d)(\\d)", "\\2", all_phis$Param
  #   )
  # )
  # phi_mat <- matrix(0, ncol = infos$q * infos$maxLag, nrow = infos$q)

  # loop across all phis
  # n_phi <- nrow(all_phis)
  phi_mat_vec <- c()
  # this loop is wrinkeling my brain
  for (i in 1:infos$q) {
    for (j in 1:infos$q) {
      for (k in 1:infos$maxLag) {
        # if phi(j)_ik exists in model frame (i.e., is not fixed to zero),
        # paste in vector
        if (any(all_phis$Param == paste0("phi(", k, ")_", i, j))) {
          if (j == infos$q & k == infos$maxLag) {
            # linebreak after last cell in row
            phi_mat_vec <- c(
              phi_mat_vec, paste0(
                "\\phi_{(", k, ")", i, j, ifelse(
                  all_phis[
                    all_phis$Param == paste0("phi(", k, ")_", i, j), "isRandom"
                  ] == 1,
                  ",i", ""
                ),
                "} \\\\"
              )
            )
          } else {
            phi_mat_vec <- c(
              phi_mat_vec, paste0(
                "\\phi_{(", k, ")", i, j, ifelse(
                  all_phis[
                    all_phis$Param == paste0("phi(", k, ")_", i, j), "isRandom"
                  ] == 1,
                  ",i", ""
                ),
                "} & "
              )
            )
          }
        } else {
          # if phi(j)_ik does not exist in model frame, paste empty cell
          if (j == infos$q & k == infos$maxLag) {
            # linebreak after last cell in row
            phi_mat_vec <- c(
              phi_mat_vec, "0 \\\\"
            )
          } else {
            phi_mat_vec <- c(
              phi_mat_vec, "0 & "
            )
          }
        }
      }
    }
  }
  # paste together
  phi_mat <- paste(phi_mat_vec, collapse = "\n")

  # initiate empty time series vector
  ts_vec <- c()
  # time series loop
  for (i in 1:infos$q) {
    for (j in 1:infos$maxLag) {
      ts_vec <- c(
        ts_vec, paste0(
          ifelse(infos$isLatent == TRUE, "\\eta_{", "y_{"),
          i, ",i(t - ", j, ")}^w \\\\"
        )
      )
    }
  }
  ts <- paste(ts_vec, collapse = "\n")

  # initiate empty innovation vector
  innos_vec <- c()
  # innovations loop
  for (i in 1:infos$q) {
    innos_vec <- c(innos_vec, paste0("\\zeta_{", i, ", it} \\\\"))
  }
  innos <- paste(innos_vec, collapse = "\n")

  # if innovation covariances are random: add eta factor to model
  if (any(model[startsWith(model$Param, "ln.sigma_"), "isRandom"] == 1)) {
    # eta loading matrix
    eta_load <- paste(infos$inno_cov_load, collapse = " \\\\\n")
    # eta matrix
    eta_inno <- "\\eta_{\\zeta_{12}, it}"
  }

  # create right hand side formula
  wmf_rhs <- paste(
    begin_bmatrix, phi_mat, end_bmatrix,
    begin_bmatrix, ts, end_bmatrix, "+",
    if (any(model[startsWith(model$Param, "ln.sigma_"), "isRandom"] == 1)) {
      paste(
        begin_bmatrix, eta_load, end_bmatrix,
        begin_bmatrix, eta_inno, end_bmatrix, "+",
        collapse = "\n"
      )
    },
    begin_bmatrix, innos, end_bmatrix,
    collapse = "\n"
  )

  if (infos$q == 1) {
    inno_dist <- paste0(
      ",\\\\ \n\\text{with}~",
      "\\zeta_{1,it} \\sim \\mathit{N}(0, \\sigma^2_{\\zeta_{1}", ifelse(
        model[model$Param == "ln.sigma2_1", "isRandom"] == 1,
        ",i", ""
      ),"})"
    )
    psi_mat_vec <- NULL
  } else if (infos$q > 1 & infos$n_inno_cors == 0 & infos$n_inno_covs == 0) {
    inno_dist <- paste0(
      ",\\\\ \n\\text{with}~",
      "\\zeta_{k,it} \\sim \\mathit{N}(\\mathbf{0}, \\sigma^2_{\\zeta_{k}",
      ifelse(
        # add i index if there are random elements in sigmas
        any(model[
          grepl("Fix", model$Type) & grepl("sigma", model$Param), "isRandom"
        ] == 1),
        ",i", ""
      ), "})"
    )
    psi_mat_vec <- NULL
  } else if (infos$n_inno_covs > 0) {
    # in case of random innovation covariance, include eta distribution
    # together with innovation distribution
    inno_dist <- paste0(
      ", \\\\ \n\\text{with}~",
      "\\zeta_{k,it} \\sim \\mathit{N}(\\mathbf{0}, \\sigma^2_{\\zeta_{k}",
      ifelse(
        # add i index if there are random elements in sigmas
        any(model[
          grepl("Fix", model$Type) & grepl("sigma", model$Param), "isRandom"
        ] == 1),
        ",i", ""
      ), "})",
      ", \\\\ \n\\text{and}~\\eta_{\\zeta_{12},it} \\sim \\mathit{N}(",
      "0, \\sigma_{\\zeta_{12}, i})"
    )
    psi_mat_vec <- NULL
  } else {
    inno_dist <- paste0(
      ",\\\\ \n\\text{with}~",
      "\\zeta_{it} \\sim \\mathit{MVN}(\\mathbf{0}, \\mathbf{\\Psi}",
      ifelse(
        # add i index if there are random elements in PSI
        any(model[
          grepl("Fix", model$Type) & grepl("sigma", model$Param), "isRandom"
        ] == 1),
        "_i", ""
      ), ")"
    )
      # psi matrix
    psi_mat_vec <- c()
    for (i in 1:infos$q) {
      for (j in 1:infos$q) {
        if (i == j) {
          if (j == infos$q) {
            psi_mat_vec <- c(psi_mat_vec, paste0(
              "\\sigma^2_{\\zeta_{", i, "}", ifelse(
                model[model$Param == paste0("ln.sigma2_", i), "isRandom"] == 1,
                ",i", ""
              ), "} \\\\"
            ))
          } else {
            psi_mat_vec <- c(psi_mat_vec, paste0(
              "\\sigma^2_{\\zeta_{", i, "}", ifelse(
                model[model$Param == paste0("ln.sigma2_", i), "isRandom"] == 1,
                ",i", ""
              ), "} &"
            ))
          }
        } else if (i > j) {
          psi_mat_vec <- c(psi_mat_vec, paste0(" & "))
        } else {
          if (j == infos$q) {
            psi_mat_vec <- c(psi_mat_vec, paste0(
              "\\sigma_{\\zeta_{", i, j, "}", ifelse(
                model[model$Param == paste0("ln.sigma_", i, j), "isRandom"] == 1,
                ",i", ""
              ), "} \\\\"
            ))
          } else {
            psi_mat_vec <- c(psi_mat_vec, paste0(
              "\\sigma_{\\zeta_{", i, j, "}", ifelse(
                model[model$Param == paste0("ln.sigma_", i, j), "isRandom"] == 1,
                ",i", ""
              ), "} &"
            ))
          }
        }
      }
    }

    # paste together
    psi_mat <- paste(psi_mat_vec, collapse = "\n")
    PSI <- paste(
      ", \\\\ \n\\text{and}~\\mathbf{\\Psi}", ifelse(
        # add i index if there are random elements in PSI
        any(model[
          grepl("Fix", model$Type) & grepl("sigma", model$Param), "isRandom"
        ] == 1),
        "_i", ""
      ), " = ",
      paste(begin_bmatrix, psi_mat, end_bmatrix, collapse = "\n"),
      collapse = "\n"
    )
  }

  # within-model formula
  wmf <- paste(
    begin_math,
    wmf_lhs, "=", wmf_rhs,
    inno_dist,
    if (!is.null(psi_mat_vec)) {
      paste(PSI)
    },
    end_math
  )

  # with caption
  within_model <- paste(begin_center, wmf_caption, end_center, wmf)


  # paste within_model to markdown
  cat(within_model, file = rmd_file, append = TRUE)


  # between-model #############################################################

  # caption
  bmf_caption <- "Between-level model."

  all_bpars <- model[grepl("Fix", model$Type), ]
  # replace dots with nothing
  all_bpars$Param <- gsub("\\.", "", all_bpars$Param)
  # replace rzeta with psi (ugly fix but ok)
  all_bpars$Param <- gsub("rzeta", "sigma", all_bpars$Param)
  all_bpars$names <- gsub(
    # replace underscore digit with latex subscript
    "(\\w+)_(\\d+)", "\\\\\\1_{\\2}", gsub(
      # replace lnsigma with ln(sigma_zeta)
      "(lnsigma2)_(\\d+)", "\\\\ln(\\\\sigma^2_{\\\\zeta_{\\2}})", gsub(
        # replace uppercase B for between-level latent variables
        "(\\w+)(B)_(\\d)", "\\\\\\1^{b}_{\\3}", gsub(
          # replace underscore digit with latex subscript for phis
          "(\\w+)(\\(\\d\\))_(\\d+)", "\\\\\\1_{\\2\\3}", gsub(
            # replace rzeta with psi in case of fixed
            # innovation covariance (ugly fix but ok)
            "(sigma)_(\\d+)", "\\\\sigma_{\\\\zeta_{\\2}}", gsub(
              # replace lnsigma12 (random innovation covariance)
              # with ln(psi)
              "(lnsigma)_(\\d+)", "\\\\ln(\\\\sigma_{\\\\zeta_{\\2}})", all_bpars$Param
            )
          )
        )
      )
    )
  )
  # special case for fixed innovation variances
  all_bpars$names <- ifelse(
    startsWith(all_bpars$Param, "sigma_"),
    gsub("(sigma)_(\\d+)", "\\\\sigma^2_{\\\\zeta_{\\2}}", all_bpars$Param),
    all_bpars$names
  )
  # add one unit-specific names version (add i subscript)
  all_bpars$names2 <- gsub(
    # replace underscore digit with latex subscript
    "(\\w+)_(\\d+)", "\\\\\\1_{\\2,i}", gsub(
      # replace lnsigma with ln(sigma_zeta)
      "(lnsigma2)_(\\d+)", "\\\\ln(\\\\sigma^2_{\\\\zeta_{\\2},i})", gsub(
        # replace uppercase B for between-level latent variables
        "(\\w+)(B)_(\\d)", "\\\\\\1^{b}_{\\3,i}", gsub(
          # replace underscore digit with latex subscript for phis
          "(\\w+)(\\(\\d\\))_(\\d+)", "\\\\\\1_{\\2\\3,i}", gsub(
            # replace rzeta with psi in case of fixed
            # innovation covariance (ugly fix but ok)
            "(sigma)_(\\d+)", "\\\\sigma_{\\\\zeta_{\\2}, i}", gsub(
              # replace lnsigma12 (random innovation covariance)
              # with ln(psi)
              "(lnsigma)_(\\d+)", "\\\\ln(\\\\sigma_{\\\\zeta_{\\2}, i})", all_bpars$Param
            )
          )
        )
      )
    )
  )
  # special case for fixed innovation variances
  all_bpars$names2 <- ifelse(
    startsWith(all_bpars$Param, "sigma_"),
    gsub("(sigma)_(\\d+)", "\\\\sigma_{\\\\zeta_{\\2},i}", all_bpars$Param),
    all_bpars$names2
  )
  # number of bpars
  n_bpars <- nrow(all_bpars)
  # left-hand side
  bpars_vec <- c()
  for (i in 1:n_bpars) {
    bpars_vec <- c(bpars_vec, paste0(all_bpars$names2[i], "\\\\"))
  }
  bpars <- paste(bpars_vec, collapse = "\n")
  # create between model left hand side
  bmf_lhs <- paste(begin_bmatrix, bpars, end_bmatrix, collapse = "\n")

  # right-hand side fixed effects
  fixef_vec <- c()
  for (i in 1:n_bpars) {
    fixef_vec <- c(fixef_vec, paste0("\\gamma_{0,", all_bpars$names[i], "}\\\\"))
  }
  fixef <- paste(fixef_vec, collapse = "\n")

  # right-hand side random effects
  ranef_vec <- c()
  for (i in 1:n_bpars) {
    if (all_bpars$isRandom[i] == 1) {
      ranef_vec <- c(ranef_vec, paste0("\\upsilon_{", all_bpars$names[i], ",i}\\\\"))
    } else {
      ranef_vec <- c(ranef_vec, paste0("0\\\\"))
    }
  }
  ranef <- paste(ranef_vec, collapse = "\n")

  # for predicted random effects
  if (nrow(infos$RE.PREDS) > 0) {
    # create vector of predictors
    re_pred_vars_vec <- c()
    n_re_preds <- length(infos$n_cov_vars)
    for (i in 1:n_re_preds) {
      re_pred_vars_vec <- c(re_pred_vars_vec, paste0(infos$n_cov_vars[i], "\\\\"))
    }
    re_pred_vars <- paste(re_pred_vars_vec, collapse = "\n")

    # create coefficient matrix
    re_preds <- infos$RE.PREDS[, c(
      "re_as_dv", "re_preds", "pred_no"
    )]
    # rename for easier binding later
    names(re_preds)[1] <- "dv"
    names(re_preds)[2] <- "pred"
    # rename parameters in dv column
    re_preds$dv <- gsub("\\.", "", re_preds$dv)
    re_preds$names <- gsub(
      # replace underscore digit with latex subscript
      "(\\w+)_(\\d+)", "\\\\\\1_{\\2}", gsub(
        # replace lnsigma with ln(pi)
        "(lnsigma2)_(\\d+)", "\\\\ln(\\\\sigma^2_{\\\\zeta_{\\2}})", gsub(
          # replace uppercase B for between-level latent variables
          "(\\w+)(B)_(\\d)", "\\\\\\1^{b}_\\3", gsub(
            # replace underscore digit with latex subscript for phis
            "(\\w+)(\\(\\d\\))_(\\d+)", "\\\\\\1_{\\2\\3}", gsub(
              # replace rzeta with psi in case of fixed
              # innovation covariance (ugly fix but ok)
              "(sigma)_(\\d+)", "\\\\sigma_{\\\\zeta_{\\2}}", gsub(
                # replace lnsigma12 (random innovation covariance)
                # with ln(psi)
                "(lnsigma)_(\\d+)", "\\\\ln(\\\\sigma_{\\\\zeta_{\\2}})", re_preds$dv
              )
            )
          )
        )
      )
    )
    # pivot to wide to create predictor matrix
    re_preds_wide <- stats::reshape(
      data = re_preds,
      idvar = "dv",
      direction = "wide",
      timevar = "pred_no",
      v.names = c("names", "pred")
    )
    # create coefficient matrix and vector
    re_coefs_mat <- matrix(0, nrow = n_bpars, ncol = n_re_preds)
    re_coefs_vec <- c()
    for (i in 1:n_bpars) {
      for (j in 1:n_re_preds) {
        re_coefs_mat[i, j] <- ifelse(
          !is.na(re_preds_wide[i, paste0("names.", j)]),
          paste0(
            "\\gamma_{", j, ",", re_preds_wide[i, paste0("names.", j)], "}"
          ),
          # paste a 0 if parameter is fixed to 0
          "0"
        )
      }
      re_coefs_vec[i] <- paste(re_coefs_mat[i, ], collapse = " & ")
    }
    # paste together
    re_coefs <- paste(re_coefs_vec, collapse = "\\\\\n")

    # paste in tex
    re_predictors <- paste(
      begin_bmatrix, re_coefs, end_bmatrix,
      begin_bmatrix, re_pred_vars, end_bmatrix,
      "+",
      collapse = "\n"
    )
  }

  # for predicted outcomes
  if (nrow(infos$OUT) > 0) {
    # create vector of outcomes and residuals
    out_vars_vec <- c()
    out_res_vec <- c()
    n_out <- infos$n_out
    for (i in 1:n_out) {
      out_vars_vec <- c(out_vars_vec, paste0(infos$out_var[i], "\\\\"))
      out_res_vec <- c(out_res_vec, paste0("\\varepsilon_{", infos$out_var[i], "}\\\\"))
    }
    out_vars <- paste(out_vars_vec, collapse = "\n")
    out_res <- paste(out_res_vec, collapse = "\n")

    # create coefficient matrix
    out <- infos$OUT[, c(
      "Var", "Pred", "out_var_no"
    )]
    # rename for easier binding later
    names(out)[1] <- "dv"
    names(out)[2] <- "pred"
    # rename parameters in dv column
    out$pred <- gsub("\\.", "", out$pred)
    out$names <- gsub(
      # replace underscore digit with latex subscript
      "(\\w+)_(\\d+)", "\\\\\\1_{\\2}", gsub(
        # replace lnsigma with ln(sigma_zeta)
        "(lnsigma2)_(\\d+)", "\\\\ln(\\\\sigma^2_{\\\\zeta_{\\2}})", gsub(
          # replace uppercase B for between-level latent variables
          "(\\w+)(B)_(\\d)", "\\\\\\1^{b}_\\3", gsub(
            # replace underscore digit with latex subscript for phis
            "(\\w+)(\\(\\d\\))_(\\d+)", "\\\\\\1_{\\2\\3}", gsub(
              # replace rzeta with psi in case of fixed
              # innovation covariance (ugly fix but ok)
              "(sigma)_(\\d+)", "\\\\sigma_{\\\\zeta_{\\2}}", gsub(
                # replace lnsigma12 (random innovation covariance)
                # with ln(psi)
                "(lnsigma)_(\\d+)", "\\\\ln(\\\\sigma_{\\\\zeta_{\\2}})", out$pred
              )
            )
          )
        )
      )
    )
    # pivot to wide to create predictor matrix
    out_wide <- stats::reshape(
      data = out,
      idvar = "pred",
      direction = "wide",
      timevar = "out_var_no",
      v.names = c("names", "dv")
    )

    # create vector of predictors (= random effects)
    out_preds <- paste(unique(out$names), collapse = "\\\\\n")
    n_out_preds <- nrow(out_wide)

    # create coefficient matrix and vector
    out_coefs_mat <- matrix(0, nrow = n_out, ncol = n_out_preds)
    out_coefs_vec <- c()
    for (i in 1:n_out) {
      for (j in 1:n_out_preds) {
        out_coefs_mat[i, j] <- ifelse(
          !is.na(out_wide[j, paste0("dv.", i)]),
          paste0(
            "\\gamma_{", j, ",", out_wide[j, paste0("dv.", i)], "}"
          ),
          # paste a 0 if parameter is fixed to 0
          "0"
        )
      }
      out_coefs_vec[i] <- paste(out_coefs_mat[i, ], collapse = " & ")
    }

    # paste together
    out_coefs <- paste(out_coefs_vec, collapse = "\\\\\n")

    # paste in tex
    out_predictors <- paste(
      begin_bmatrix, out_vars, end_bmatrix, "=",
      begin_bmatrix, out_coefs, end_bmatrix,
      begin_bmatrix, out_preds, end_bmatrix, "+",
      begin_bmatrix, out_res, end_bmatrix,
      collapse = "\n"
    )
  }

  # create right hand side formula
  bmf_rhs <- paste(
    begin_bmatrix, fixef, end_bmatrix, "+\n",
    # paste predictors here if provided
    if (nrow(infos$RE.PREDS) > 0) {paste(re_predictors)},
    begin_bmatrix, ranef, end_bmatrix,
    collapse = "\n"
  )

  ranef_dist <- ",\\\\ \n\\text{with}~
  \\upsilon_{i} \\sim \\mathit{MVN}(\\mathbf{0}, \\mathbf{\\Omega})"

  # between-model formula
  bmf <- paste(
    begin_math,
    bmf_lhs, "=", bmf_rhs,
    ranef_dist,
    # paste outcome prediction here if provided
    if (nrow(infos$OUT) > 0) {paste("\\\\\n", out_predictors)},
    end_math
  )

  # with caption
  between_model <- paste(begin_center, bmf_caption, end_center, bmf)


  # paste within_model to markdown
  cat(between_model, file = rmd_file, append = TRUE)

  # render markdown input #####################################################

  # render markdown
  rmarkdown::render(
    input = rmd_file
  )

}
