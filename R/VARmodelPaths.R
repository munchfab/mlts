#' Create Path Diagrams from VARmodel object
#'
#' @param VARmodel The VARmodel object.
#' @param data An optional `data.frame` including the variables to be used in
#' the path diagram.
#' @param labels An optional character string inluding the names to be used in
#' the path diagram.
#' @param add.png Logical. Set to `TRUE` to transform created PDF as .png files
#'  using `pdftools::pdf_convert`.
#' @return An Rmarkdown file that is automatically rendered to a pdf document.
#' @export
#'
#' @examples
VARmodelPaths <- function(VARmodel, data = NULL, labels = NULL, add.png = FALSE) {


  # extrat model infos
  infos <- VARmodelEval(VARmodel)

  # extract model data frame
  model <- VARmodel

  # create empty markdown file, delete if already existing
  if (file.exists("pathmodel.rmd")) {
    file.remove("pathmodel.rmd")
  }
  rmarkdown::draft(file = "pathmodel.rmd",
                   template = "pathmodel",
                   package = "dsemr",
                   edit = FALSE)

  # latex figure begin and end
  begin_figure <- "\n\\begin{figure}\n\\centering"
  end_figure <- "\\end{figure}\n"

  # create empty string with tikz specifications
  begin_tikz <- "\\begin{tikzpicture}[
    auto, > = latex, align=center,
  	latent/.style = {
  	  circle, draw, thick, inner sep = 2pt, minimum width = \\Radius
  	},
  	manifest/.style = {
  	  rectangle, draw, thick, inner sep = 0pt, minimum size = \\Radius/2
  	},
  	intercept/.style = {
  	  regular polygon,regular polygon sides = 3, draw, thick,
  	  inner sep = 0pt, minimum size = \\Radius
  	},
  	mean/.style = {
  	  regular polygon, regular polygon sides = 3, draw, thick,
  	  inner sep = 0pt, minimum size = 8mm
  	},
  	path/.style = {
  	  arrows = ->, thick, > = {stealth[]}
  	},
  	error/.style = {
  	  circle, draw = none, fill = none, thick,
  	  inner sep = 0pt, minimum size = 5mm
  	},
  	var/.style = {
  	  <->, thick, > = {stealth[]}, bend right = 270, looseness = 2
  	},
  	cov/.style = {
  	  <->, thick, > = {stealth[]}, bend right = 300
  	},
  	% style to add a circle in the middle of a path
    random/.style = {
      postaction = {
        decorate, decoration = {
          markings,
          mark = at position .5 with {
            \\draw[fill = black] circle[radius = 2pt];
          }
        }
      }
    },
    random_cl/.style = {
      postaction = {decorate, decoration = {
        markings,
        mark = at position .2 with {
          \\draw[fill = black] circle[radius = 2pt];}
        }
      }
    },
  	]"

  # end tikz picture
  end_tikz <- "\\end{tikzpicture}\n"


  # decomposition #############################################################

  dc_caption <- "% caption
  \\node  [above = 1em, align = center]  at  (current bounding box.north)  {Decomposition.};"

  if (infos$q == 1) {
    dc <- "
    % draw decomposition
    \\node  [manifest]  (y1t)  {$y_{1,t}$};
    \\node  [latent]  (y1wt)  [above = 2.5em of y1t]  {$y_{1,t}^w$};
    \\node  [latent]  (mu1)  [below = 2.5em of y1t]  {$\\mu_{1}$};

    % draw paths
    \\draw  [path]  (y1wt)  to node  []  {}  (y1t);
    \\draw  [path]  (mu1)  to node  []  {}  (y1t);
    "
    if (infos$p > 1) {
      # initiate empty vectors to store tikz nodes
      ind_vec <- c() # indicator variables
      # latent variables (within and between)
      wlat <- "\\node  [latent]  (eta1wt)  [above = 3.5em of c]  {$\\eta^w_{1,t}$};"
      blat <- "\\node  [latent]  (mu1)  [below = 3.5em of c]  {$\\mu_{1}$};"
      wlat_paths_vec <- c() # paths from latent variables to indicators
      blat_paths_vec <- c() # paths from latent variables to indicators
      eps_vec <- c() # residuals
      eps_paths_vec <- c() # paths from residuals to indicators
      for (i in 1:infos$p) {
        eps_vec <- c(eps_vec, paste0(
          "\\node  [error]  (eps", i, "t)  [above left = .75em of y", i, "t]  {};"
        ))
        if (i == 1) {
          # first indicator variable
          ind_vec <- c(ind_vec, paste0(
            "\\node  [manifest]  (y", i, "t)  {$y_{", i, ",t}$};"
          ))
          # first error (labeled)
          eps_vec <- c(eps_vec, paste0(
            "\\node  [error]  (eps", i, "t)  [above left = 1em of y",
            i, "t]  {\\scriptsize$\\varepsilon^w_{", i, ",t}$};"
          ))
        } else {
          # all other indicators
          ind_vec <- c(ind_vec, paste0(
            "\\node  [manifest]  (y", i, "t)  [right = 1.5em of y",
            i - 1, "t]  {$y_{", i, ",t}$};"
          ))
          if (i == infos$p) {
            # place coordinates between first and last indicator variables
            # of construct
            coord <- paste0(
              "\\node  [coordinate]  (c)  at ($(y1t) !0.5! (y", i, "t)$)  {};"
            )
          }
        }
        # draw paths
        wlat_paths_vec <- c(wlat_paths_vec, paste0(
          "\\draw  [path]  (eta1wt)  to node  [fill = white, anchor = center]
            {\\scriptsize$\\lambda^w_{", i, "}$}  (y", i, "t.north);"
        ))
        blat_paths_vec <- c(blat_paths_vec, paste0(
          "\\draw  [path]  (mu1)  to node  [fill = white, anchor = center]
            {\\scriptsize$\\lambda^{\\mu}_{", i, "}$}  (y", i, "t.south);"
        ))
        eps_paths_vec <- c(eps_paths_vec, paste0(
          "\\draw  [path]  (eps", i, "t)  to node  []  {}  (y", i, "t.north west);"
        ))
      }
      ind <- paste(ind_vec, collapse = "\n")
      wlat_paths <- paste(wlat_paths_vec, collapse = "\n")
      blat_paths <- paste(blat_paths_vec, collapse = "\n")
      eps <- paste(eps_vec, collapse = "\n")
      eps_paths <- paste(eps_paths_vec, collapse = "\n")
      # paste together
      dc <- paste(ind, coord, wlat, blat,
                  wlat_paths, blat_paths,
                  eps, eps_paths, sep = "\n")
    }
  } else { # for q > 1
    dc <- paste0(
    "
    % draw decomposition
    \\node  [manifest]  (y1t)  {$y_{1,t}$};
    \\node  [latent]  (y1wt)  [above = 2.5em of y1t]  {$y_{1,t}^w$};
    \\node  [latent]  (mu1)  [below = 2.5em of y1t]  {$\\mu_{1}$};

    % draw nodes
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    infos$q, "}
    \\node  [manifest]  (y\\i t)  [right = 2.5em of y\\lasti t]  {$y_{\\i,t}$};

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    infos$q, "}
    \\node  [latent]  (y\\i wt)  [above = 2.5em of y\\i t]  {$y_{\\i,t}^w$};

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ...,",
    infos$q, "}
    \\node  [latent]  (mu\\i)  [below = 2.5em of y\\i t]  {$\\mu_{\\i}$};

    % draw paths
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
    infos$q, "}
    \\draw  [path]  (y\\i wt)  to node  []  {}  (y\\i t);

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
    infos$q, "}
    \\draw  [path]  (mu\\i)  to node  []  {}  (y\\i t);"
    )
    if (any(infos$p > 1)) {
      # initiate empty vectors to store tikz nodes
      ind_vec <- c() # indicator variables
      coord_vec <- c() # coordinates for positioning of latent variables
      wlat_vec <- c() # within-level latent variables
      blat_vec <- c() # between-level latent variables
      wlat_paths_vec <- c() # paths from latent variables to indicators
      blat_paths_vec <- c() # paths from latent variables to indicators
      eps_vec <- c() # residuals
      eps_paths_vec <- c() # paths from residuals to indicators
      for (i in 1:infos$q) {
        wlat_vec <- c(wlat_vec, paste0(
          "\\node  [latent]  (etaw", i, "t)  [above = 3.5em of c", i,
          "]  {$\\eta^w_{", i, ",t}$};"
        ))
        blat_vec <- c(blat_vec, paste0(
          "\\node  [latent]  (etab", i, ")  [below = 3.5em of c", i,
          "]  {$\\mu_{", i, "}$};"
        ))
        for (j in 1:infos$p[i]) {
          eps_vec <- c(eps_vec, paste0(
            "\\node  [error]  (eps", i, j, "t)  [above left = .75em of y", i, j, "t]  {};"
          ))
          if (i == 1 & j == 1) {
            # first indicator variable
            ind_vec <- c(ind_vec, paste0(
              "\\node  [manifest]  (y", i, j, "t)  {$y_{", i, j, ",t}$};"
            ))
            # first error (labeled)
            eps_vec <- c(eps_vec, paste0(
              "\\node  [error]  (eps", i, j, "t)  [above left = 1em of y",
              i, j, "t]  {\\scriptsize$\\varepsilon^w_{", i, j, ",t}$};"
            ))
          } else if (i > 1 & j == 1) {
            # first indicator variable of next construct
            ind_vec <- c(ind_vec, paste0(
              "\\node  [manifest]  (y", i, j, "t)  [right = 1.5em of y",
              i - 1, infos$p[i - 1], "t]  {$y_{", i, j, ",t}$};"
            ))
          } else {
            # all other indicators
            ind_vec <- c(ind_vec, paste0(
              "\\node  [manifest]  (y", i, j, "t)  [right = 1.5em of y",
              i, j - 1, "t]  {$y_{", i, j, ",t}$};"
            ))
            if (j == infos$p[i]) {
              # place coordinates between first and last indicator variables
              # of construct i
              coord_vec <- c(coord_vec, paste0(
                "\\node  [coordinate]  (c", i,
                ")  at ($(y", i, "1t) !0.5! (y", i, j, "t)$)  {};"
              ))
            }
          }
          # draw paths
          wlat_paths_vec <- c(wlat_paths_vec, paste0(
            "\\draw  [path]  (etaw",
            i, "t)  to node  [fill = white, anchor = center]  {\\scriptsize$\\lambda^w_{",
            i, j, "}$}  (y", i, j, "t.north);"
          ))
          blat_paths_vec <- c(blat_paths_vec, paste0(
            "\\draw  [path]  (etab",
            i, ")  to node  [fill = white, anchor = center]  {\\scriptsize$\\lambda^{\\mu}_{",
            i, j, "}$}  (y", i, j, "t.south);"
          ))
          eps_paths_vec <- c(eps_paths_vec, paste0(
            "\\draw  [path]  (eps", i, j, "t)  to node  []  {}  (y", i, j, "t.north west);"
          ))
        }
      }
      ind <- paste(ind_vec, collapse = "\n")
      coord <- paste(coord_vec, collapse = "\n")
      wlat <- paste(wlat_vec, collapse = "\n")
      blat <- paste(blat_vec, collapse = "\n")
      wlat_paths <- paste(wlat_paths_vec, collapse = "\n")
      blat_paths <- paste(blat_paths_vec, collapse = "\n")
      eps <- paste(eps_vec, collapse = "\n")
      eps_paths <- paste(eps_paths_vec, collapse = "\n")
      # paste together
      dc <- paste(ind, coord, wlat, blat,
                  wlat_paths, blat_paths,
                  eps, eps_paths, sep = "\n")
    }
  }

  # paste together
  decomposition <- paste(
    # begin_figure, # not needed
    # dc_caption,
    begin_tikz,
    dc,
    dc_caption,
    end_tikz,
    # end_figure,
    sep = "\n"
  )

  # paste decomposition to markdown
  cat(decomposition, file = "pathmodel.rmd", append = TRUE)


  # within-model ##############################################################

  # within-model caption
  wm_caption <- "% caption
  \\node  [above = 1em, align = center]  at  (current bounding box.north)  {Within-model.};"

  # draw nodes for q time-series constructs
  if (infos$q == 1) {
    wm <- "
    % draw within-level structural model
    \\node  [latent]  (y1wt-1)  {$y_{1,t-1}^w$};
    \\node  [latent]  (y1wt)  [right = 5em of y1wt-1]  {$y_{1,t}^w$};
    \\node  [latent]  (delta1t)  [right = 2.5em of ywt]  {$\\delta_{{y_1},t}$};

    % draw paths
    \\draw  [path]  (deltat)  to node  []  {}  (ywt);
    "
    if (infos$p > 1) {
      wm <- "
      % draw within-level structural model
      \\node  [latent]  (eta1wt-1)  {$\\eta_{1,t-1}^w$};
      \\node  [latent]  (eta1wt)  [right = 5em of eta1wt-1]  {$\\eta_{1,t}^w$};
      \\node  [latent]  (delta1t)  [right = 2.5em of eta1wt]  {$\\delta_{{\\eta^w_{1}},t}$};

      % draw paths
      \\draw  [path]  (delta1t)  to node  []  {}  (eta1wt);
      "
    }
    # draw paths conditional on isRandom
    phi <- paste0(
      "\\draw  [path", ifelse(
        model[
          model$Param == "phi_11" & model$Type == "Fix effect", "isRandom"
        ] == 1,
        ", postaction = random]", "]"
      ), ifelse(
        # select start node conditional on infos$p
        infos$p > 1, "  (eta1wt-1)", "  (y1wt-1)"
      ), "  to node  []  {$\\phi_{11}$}  (", ifelse(
        # select target node conditional on infos$p
        infos$p > 1, "eta1wt", "y1wt"
      ), ");"
    )
    ln.sigma2 <- paste0(
      "\\draw  [var", ifelse(
        model[
          model$Param == "ln.sigma2_1" & model$Type == "Fix effect", "isRandom"
        ] == 1,
        ", postaction = random]", "]"
      ),
      "  (delta1t.120)  to node  []  {$\\pi_{1}$}  (delta1t.60);"
    )
  } else { # for q > 1
    wm <- paste0(
    "
    % draw within-level structural model
    \\node  [latent]  (y1wt-1)  {$y_{1,t-1}^w$};
    \\node  [latent]  (y1wt)  [right =5em of y1wt-1]  {$y_{1,t}^w$};
    \\node  [latent]  (delta1t)  [right =2.5em of y1wt]  {$\\delta_{1,t}$};

    % draw nodes
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
      infos$q, "}
    \\node  [latent]  (y\\i wt-1)  [below = 2.5em of y\\lasti wt-1]  {$y_{\\i ,t-1}^w$};

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    infos$q, "}
    \\node  [latent]  (y\\i wt)  [below = 2.5em of y\\lasti wt]  {$y_{\\i ,t}^w$};

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    infos$q, "}
    \\node  [latent]  (delta\\i t)  [below = 2.5em of delta\\lasti t]  {$\\delta_{\\i ,t}$};

    % draw paths from residuals to yiwt
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
    infos$q, "}
    \\draw  [path]  (delta\\i t)  to node  []  {}  (y\\i wt);
    ")
    if (any(infos$p > 1)) {
      wm <- paste0(
      "
      % draw within-level structural model
      \\node  [latent]  (eta1wt-1)  {$\\eta_{1,t-1}^w$};
      \\node  [latent]  (eta1wt)  [right = 5em of eta1wt-1]  {$\\eta_{1,t}^w$};
      \\node  [latent]  (delta1t)  [right = 2.5em of eta1wt]  {$\\delta_{1,t}$};

      % draw nodes
      \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
          infos$q, "}
      \\node  [latent]  (eta\\i wt-1)  [below = 2.5em of eta\\lasti wt-1]  {$\\eta_{\\i ,t-1}^w$};

      \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
          infos$q, "}
      \\node  [latent]  (eta\\i wt)  [below = 2.5em of eta\\lasti wt]  {$\\eta_{\\i ,t}^w$};

      \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
          infos$q, "}
      \\node  [latent]  (delta\\i t)  [below = 2.5em of delta\\lasti t]  {$\\delta_{\\i ,t}$};

      % draw paths from residuals to etaiwt
      \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
          infos$q, "}
      \\draw  [path]  (delta\\i t)  to node  []  {}  (eta\\i wt);
      ")
    }

    # store number of phi-parameters for loops
    all_phis <- model[grepl("phi", model$Param) & grepl("Fix", model$Type), ]
    n_phi <- nrow(all_phis)
    phi_vec <- c()
    # draw paths conditional on isRandom
    # and do that for every phi in the model
    for (i in 1:sqrt(n_phi)) {
      # loop across phi subscripts
      for (j in 1:sqrt(n_phi)){
        # for autoregressive paths
        if (i == j) {
          phi_vec <- c(
            phi_vec, paste0(
              "\\draw  [path", ifelse(
                # decorate with dot on path if parameter is random
                all_phis[all_phis$Param == paste0("phi_", i, j), "isRandom"] == 1,
                ", postaction = random]", "]"
              ),
              "  (", ifelse(
                # select start node conditional on infos$p
                infos$p > 1, "eta", "y"
              ), i, "wt-1)  to node  []  {$\\phi_{", i, j, "}$}  (", ifelse(
                # select target node conditional on infos$p
                infos$p > 1, "eta", "y"
              ), i, "wt);\n"
            )
          )
        } else {
          # for cross-lagged paths
          phi_vec <- c(
            phi_vec, paste0(
              "\\draw  [path", ifelse(
                # decorate with dot on path if parameter is random
                all_phis[all_phis$Param == paste0("phi_", i, j), "isRandom"] == 1,
                ", postaction = random_cl]", "]"
              ),
              "  (", ifelse(
                # select start node conditional on infos$p
                infos$p > 1, "eta", "y"
              ), i, "wt-1)  to node  [pos = .2]  {$\\phi_{", i, j, "}$}  (", ifelse(
                # select target node conditional on infos$p
                infos$p > 1, "eta", "y"
              ), j, "wt);\n"
            )
          )
        }
      }
    }
    # paste phis in one string
    phi <- paste(phi_vec, collapse = "")

    # do the same for innovation variances
    all_sigmas <- model[grepl("ln.sigma2", model$Param) & grepl("Fix", model$Type), ]
    n_sigma <- nrow(all_sigmas)
    sigma_vec <- c()
    for (i in 1:n_sigma) {
      sigma_vec <- c(
        sigma_vec, paste0(
          "\\draw  [var", ifelse(
            # decorate with dot on path if parameter is random
            all_sigmas[all_sigmas$Param == paste0("ln.sigma2_", i), "isRandom"] == 1,
            ", postaction = random]", "]"
          ),
          "  (delta", i, "t.120)  to node  []  {$\\pi_{", i, "}$}  (delta", i, "t.60);\n"
        )
      )
    }
    # do the same for innovation covariances
    all_psis <- model[grepl("ln.sigma_", model$Param) & grepl("Fix", model$Type), ]
    # n_psi <- nrow(all_psis)
    psi_vec <- c()
    if (infos$q == 1) {
      psi_vec <- paste0(
        "\\draw  [cov", ifelse(
          # decorate with dot on path if parameter is random
          all_psis[all_psis$Param == "ln.sigma_12", "isRandom"] == 1,
          ", postaction = random]", "]"
        ),
        "  (delta1t.0)  to node  []  {$\\pi_{12}$}  (delta2t.0);\n"
      )
    } else {
      for (i in 1:(infos$q - 1)) {
        for (j in (i + 1):infos$q) {
          psi_vec <- c(
            psi_vec, paste0(
              "\\draw  [cov", ifelse(
                # decorate with dot on path if parameter is random
                all_psis[all_psis$Param == paste0("ln.sigma_", i, j), "isRandom"] == 1,
                ", postaction = random]", "]"
              ),
              "  (delta", i, "t.0)  to node  []  {$\\pi_{", i, j, "}$}  (delta", j, "t.0);\n"
            )
          )
        }
      }
    }
    # paste sigmas in one string
    ln.sigma2 <- paste(sigma_vec, psi_vec, collapse = "")
  }
  # paste together
  within_model <- paste(wm, phi, ln.sigma2, sep = "\n")

  # paste together
  within_model <- paste(
    # begin_figure,
    # wm_caption,
    begin_tikz,
    within_model,
    wm_caption,
    end_tikz,
    # end_figure,
    sep = "\n"
  )

  # paste within_model to markdown
  cat(within_model, file = "pathmodel.rmd", append = TRUE)


  # between-model #############################################################

  # between-model caption
  bm_caption <- "% caption
  \\node  [above = 1em, align = center]  at  (current bounding box.north)  {Between-model.};"

  # draw nodes for q time-series constructs
  if (infos$q == 1) {
    bm <- paste0("
      % draw between-level structural model
      \\node  [latent", ifelse(
        model[
          model$Param == "mu_1" & grepl("Fix", model$Type), "isRandom"
        ] == 0, ", color = gray]", "]"
      ), "  (mu_1)  {$\\mu_{1}$};
      \\node  [latent", ifelse(
        model[
          model$Param == "phi_11" & grepl("Fix", model$Type), "isRandom"
        ] == 0, ", color = gray]", "]"
      ), "  (phi_11)  [right = 1.5em of mu_1]  {$\\phi_{11}$};
      \\node  [latent", ifelse(
        model[
          model$Param == "ln.sigma2_1" & grepl("Fix", model$Type), "isRandom"
        ] == 0, ", color = gray]", "]"
      ), "  (lnsigma2_1)  [right = 1.5em of phi_11]  {$\\log(\\pi_{1})$};"
    )
    # draw (co-)variances conditional on existing correlations
    r_mu_1.phi_11 <- ifelse(
      nrow(model[model$Param == "r_mu_1.phi_11"]) > 0,
      "\\draw  [cov]  (mu_1.north)  to node  []  {}  (phi_11.north);",
      ""
    )
    r_mu_1.ln.sigma2_1 <- ifelse(
      nrow(model[model$Param == "r_mu_1.phi_11"]) > 0,
      "\\draw  [cov]  (mu_1.north)  to node  []  {}  (lnsigma2_1.north);",
      ""
    )
    r_phi_11.ln.sigma2_1 <- ifelse(
      nrow(model[model$Param == "r_mu_1.phi_11"]) > 0,
      "\\draw  [cov]  (phi_11.north)  to node  []  {}  (lnsigma2_1.north);",
      ""
    )
    # paste together
    covs <- paste(r_mu_1.phi_11, r_mu_1.ln.sigma2_1, r_phi_11.ln.sigma2_1,
                  collapse = "")

  } else { # for q > 1
    all_bpars <- model[grepl("Fix", model$Type), ]
    # the next line removes non-random effects completely
    # all_bpars <- model[grepl("Fix", model$Type) & model$isRandom == 1, ]

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
    n_bpars <- nrow(all_bpars)
    bpars_vec <- c()
    for (i in 1:n_bpars) {
      if (i == 1) {
        bpars_vec <- c(
          bpars_vec, paste0(
            "\\node  [latent", ifelse(
              # draw in gray if parameter is fixed
              all_bpars[i, "isRandom"] == 0,
              ", color = gray]  (", "]  ("
            ),
            all_bpars[i, "Param"], ")  {$\\", names_bpars[i], "$};\n"
          )
        )
      } else {
        bpars_vec <- c(
          bpars_vec, paste0(
            "\\node  [latent", ifelse(
              # draw in gray if parameter is fixed
              all_bpars[i, "isRandom"] == 0,
              ", color = gray]  (", "]  ("
            ),
            all_bpars[i, "Param"], ")  [right = 1em of ",
            all_bpars[i - 1, "Param"],
            "]  {$\\", names_bpars[i], "$};\n"
          )
        )
      }
    }
    bm <- paste(bpars_vec, collapse = "")

    # draw between-model covariances

    # place coordinates above between-level parameters
    cov_coordinates_vec <- c()
    for (i in 1:n_bpars) {
      cov_coordinates_vec <- c(
        cov_coordinates_vec, paste0(
          "\\node  [coordinate]  (c", i, ")  [above = 2em of ",
          all_bpars[i, "Param"], "]  {};\n"
        )
      )
    }
    cov_coordinates <- paste(cov_coordinates_vec, collapse = "")

    # draw arrows from coordinates to between-level parameters
    cov_arrows_vec <- c()
    for (i in 1:n_bpars) {
      # if parameter is not random, don't draw arrow
      if (all_bpars[i, "isRandom"] == 0) {
        cov_arrows_vec <- c(
          cov_arrows_vec, paste0("\n"))
      } else {
        cov_arrows_vec <- c(
          cov_arrows_vec, paste0(
            "\\draw  [path]  (c", i, ")  to node  []  {}  (",
            all_bpars[i, "Param"], ");\n"
          )
        )
      }
    }
    cov_arrows <- paste(cov_arrows_vec, collapse = "")

    # draw line above arrows to connect for covariance
    cov_line <- paste0(
      "\\draw  [thick]  (c",
      # start line above first random between-level parameter
      which(all_bpars[, "isRandom"] == 1)[1],
      ")  --  (c",
      # end line above last random between-level parameter
      which(all_bpars[, "isRandom"] == 1)[length(which(all_bpars[, "isRandom"] == 1))],
      ");\n"
    )

    covs <- paste(cov_coordinates, cov_arrows, cov_line)


    # old stuff - draw covariances separately as rounded errors

    # all_covs <- model[grepl("RE Cor", model$Param_Label), ]
    # # replace ln. with ln for easier string subsetting using dot
    # all_covs$Param <- gsub("ln\\.", "ln", gsub(
    #   # delete r_ at beginning
    #   "r_", "", all_covs$Param
    #   )
    # )
    # # this is probably not needed because paths are not labeled
    # all_covs$names <- gsub(
    #   # replace underscore digit with latex subscript
    #   "()_(\\d+)", "\\1_{\\2}", gsub(
    #     # replace lnsigma with log(pi)
    #     "(lnsigma2|lnsigma)_(\\d+)", "log(\\\\pi_{\\2})", all_covs$Param
    #   )
    # )
    # # split string at dot to generate path start and end
    # all_covs$from <- gsub("(\\w+)\\.(\\w+)", "\\1", all_covs$Param)
    # all_covs$to <- gsub("(\\w+)\\.(\\w+)", "\\2", all_covs$Param)
    # n_covs <- nrow(all_covs)
    # covs_vec <- c()
    # for (i in 1:n_covs) {
    #   covs_vec <- c(
    #     covs_vec, paste0(
    #       "\\draw  [cov]  (",
    #       all_covs$from[i],
    #       ".north)  to node  []  {}  (",
    #       all_covs$to[i], ".north);\n"
    #     )
    #   )
    # }
    # covs <- paste(covs_vec, collapse = "")
  }

  # paste together
  between_model <- paste(bm, covs, sep = "\n")

  # paste together
  between_model <- paste(
    # begin_figure,
    # bm_caption,
    begin_tikz,
    between_model,
    bm_caption,
    end_tikz,
    # end_figure,
    sep = "\n"
  )

  # paste within_model to markdown
  cat(between_model, file = "pathmodel.rmd", append = TRUE)




  # render markdown input #####################################################

  # render markdown
  rmarkdown::render(
    input = "pathmodel.rmd"
  )

  # optional: store PDF pages as separate pngs ################################
  if(add.png == T){
  pdftools::pdf_convert(
    pdf = "pathmodel.pdf",
    format = "png",
    pages = NULL,
    filenames = NULL,
    dpi = 720,
    antialias = TRUE,
    opw = "",
    upw = "",
    verbose = FALSE
  )
  }

}
