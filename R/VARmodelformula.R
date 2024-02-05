#' Create TeX Model Formula from VARmodel object
#'
#' @param VARmodel The VARmodel object.
#' @param data An optional `data.frame` including the variables to be used in the formula.
#' @param labels An optional character string inluding the names to be used in the path diagram.
#'
#' @return An Rmarkdown file that is automatically rendered to a pdf document.
#' @export
#'
#' @examples
VARmodelformula <- function(VARmodel, data = NULL, labels = NULL) {

  # extract model data frame
  model <- VARmodel$VARmodel

  # create empty markdown file, delete if already existing
  if (file.exists("formula.rmd")) {
    file.remove("formula.rmd")
  }
  rmarkdown::draft(file = "formula.rmd",
                   template = "formula",
                   package = "dsemr",
                   edit = FALSE)

  # latex center begin and end
  begin_center <- "\n\\begin{center}\n"
  end_center <- "\n\\end{center}\n"
  # latex math begin and end
  begin_math <- "\\[\n\\begin{aligned}"
  end_math <- "\n\\end{aligned}\n\\]"
  # latex bmatrix begin and end
  begin_bmatrix <- "\n\\begin{bmatrix}\n"
  end_bmatrix <- "\\end{bmatrix}"


  # within-model ##############################################################

  # caption
  wmf_caption <- "Within-model."

  # initiate empty dv vector
  dvs_vec <- c()
  # left hand side
  for (i in 1:VARmodel$q) {
    dvs_vec <- c(dvs_vec, paste0("y_{", i, ", t}^w \\\\"))
  }
  dvs <- paste(dvs_vec, collapse = "\n")
  wmf_lhs <- paste(begin_bmatrix, dvs, end_bmatrix, collapse = "\n")

  # autoregressive / cross-lagged parameter matrix

  # store number of phi-parameters for loops
  all_phis <- model[grepl("phi", model$Param) & grepl("Fix", model$Type), ]

  phi_mat_vec <- c()
  # loop across phi subscripts
  for (i in 1:VARmodel$q) {
    for (j in 1:VARmodel$q) {
      # if phi_ij exists in model frame (i.e., is not fixed to zero), paste in matrix
      if (any(all_phis$Param == paste0("phi_", i, j))) {
        if (j == VARmodel$q) {
          phi_mat_vec <- c(phi_mat_vec, paste0("\\phi_{", i, j, "} \\\\ \n"))
        } else if (i == VARmodel$q & j == VARmodel$q) {
          # line before end of bmatrix must not be broken
          phi_mat_vec <- c(phi_mat_vec, paste0("\\phi_{", i, j, "} \\\\"))
        } else {
          phi_mat_vec <- c(phi_mat_vec, paste0("\\phi_{", i, j, "} & "))
        }
      } else {
        # if phi_ij does not exist in model frame, paste empty cell
        if (j == VARmodel$q) {
          phi_mat_vec <- c(phi_mat_vec, paste0("0 \n"))
        } else if (i == VARmodel$q & j == VARmodel$q) {
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
  for (i in 1:VARmodel$q) {
    ts_vec <- c(ts_vec, paste0("y_{", i, ", t - 1}^w \\\\")
    )
  }
  ts <- paste(ts_vec, collapse = "\n")

  # initiate empty innovation vector
  innos_vec <- c()
  # innovations loop
  for (i in 1:VARmodel$q) {
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
