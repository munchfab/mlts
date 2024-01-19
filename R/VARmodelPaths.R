VARmodelPaths <- function(VARmodel, data = NULL, labels = NULL) {

  # extract model data frame
  model <- VARmodel$VARmodel

  # create empty markdown file, delete if already existing
  if (file.exists("pathmodel.rmd")) {
    file.remove("pathmodel.rmd")
  }
  rmarkdown::draft(file = "pathmodel.rmd",
                   template = "pathmodel",
                   package = "dsemr",
                   edit = FALSE)

  # latex figure begin and end
  begin_figure <- "\\begin{figure}\n\\centering"
  end_figure <- "\\end{figure}\n"

  # create empty string with tikz specifications
  begin_tikz <- "\\begin{tikzpicture}[
    auto, > = latex, align=center,
  	latent/.style = {circle, draw, thick, inner sep = 2pt, minimum width = \\Radius},
  	manifest/.style = {rectangle, draw, thick, inner sep = 0pt, minimum size = \\Radius/2},
  	intercept/.style = {regular polygon,regular polygon sides = 3, draw, thick, inner sep = 0pt, minimum size = \\Radius},
  	mean/.style = {regular polygon, regular polygon sides = 3, draw, thick, inner sep = 0pt, minimum size = 8mm},
  	path/.style = {arrows = ->, thick, > = {stealth[]}},
  	error/.style = {circle, draw = none, fill = none, thick, inner sep = 0pt, minimum size = 5mm},
  	var/.style = {<->, thick, > = {stealth[]}, bend right = 270, looseness = 2},
  	cov/.style = {<->, thick, > = {stealth[]}, bend right = 300},
  	% style to add a circle in the middle of a path
    random/.style = {postaction = {decorate, decoration = {markings, mark = at position .5 with {\\draw[fill = black] circle[radius = 2pt];}}}},
    random_cl/.style = {postaction = {decorate, decoration = {markings, mark = at position .2 with {\\draw[fill = black] circle[radius = 2pt];}}}},
  	]"

  # end tikz picture
  end_tikz <- "\\end{tikzpicture}\n"


  # decomposition #############################################################

  dc_caption <- "\\caption*{Decomposition.}"

  if (VARmodel$q == 1) {
    dc <- "
    % draw decomposition
    \\node [manifest] (y1t) {$y_{1,t}$};
    \\node [latent]   (y1wt)    [above = 2.5em of y1t]  {$y_{1,t}^w$};
    \\node [latent]   (mu_y1)   [below = 2.5em of y1t]  {$\\mu_{y_1,t}$};

    % draw paths
    \\draw [path] (y1wt) to node [] {} (y1t);
    \\draw [path] (mu_y1) to node [] {} (y1t);"
  } else {
    dc <- paste0(
    "
    % draw decomposition
    \\node [manifest] (y1t) {$y_{1,t}$};
    \\node [latent]   (y1wt)    [above = 2.5em of y1t]  {$y_{1,t}^w$};
    \\node [latent]   (mu_y1)   [below = 2.5em of y1t]  {$\\mu_{y_1,t}$};

    % draw nodes
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    VARmodel$q, "}
    \\node [manifest] (y\\i t) [right = 2.5em of y\\lasti t] {$y_{\\i,t}$};

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    VARmodel$q, "}
    \\node [latent] (y\\i wt) [above = 2.5em of y\\i t] {$y_{\\i,t}^w$};

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ...,",
    VARmodel$q, "}
    \\node [latent] (mu_y\\i) [below = 2.5em of y\\i t] {$\\mu_{y_{\\i,t}}$};

    % draw paths
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
    VARmodel$q, "}
    \\draw [path] (y\\i wt) to node [] {} (y\\i t);

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
    VARmodel$q, "}
    \\draw [path] (mu_y\\i) to node [] {} (y\\i t);
    "
    )
  }

  # paste together
  decomposition <- paste(
    begin_figure,
    dc_caption,
    begin_tikz,
    dc,
    end_tikz,
    end_figure,
    sep = "\n"
  )

  # paste decomposition to markdown
  cat(decomposition, file = "pathmodel.rmd", append = TRUE)


  # within-model ##############################################################

  # within-model caption
  wm_caption <- "\\caption*{Within-model.}"

  # draw nodes for q time-series constructs
  if (VARmodel$q == 1) {
    wm <- "
    % draw within-level structural model
    \\node [latent] (y1wt-1)  {$y_{1,t-1}^w$};
    \\node [latent] (y1wt)    [right =5em of y1wt-1]  {$y_{1,t}^w$};
    \\node [latent] (delta1t)   [right =2.5em of y1wt]  {$\\delta_{{y_1},t}$};

    % draw paths
    \\draw [path]  (delta1t)   to node [] {}           (y1wt);"
    # draw paths conditional on isRandom
    phi <- paste0(
      "\\draw [path", ifelse(
        model[model$Param == "phi_11" & model$Type == "Fix effect", "isRandom"] == 1,
        ", postaction = random]", "]"
      ),
      " (y1wt-1)  to node [] {$\\phi_{y_1}$}  (y1wt);"
    )
    ln.sigma2 <- paste0(
      "\\draw [var", ifelse(
        model[model$Param == "ln.sigma2_1" & model$Type == "Fix effect", "isRandom"] == 1,
        ", postaction = random]", "]"
      ),
      " (delta1t.120)   to node [] {$\\pi_{y_1}$} (delta1t.60);"
    )
  } else {
    wm <- paste0(
    "
    % draw within-level structural model
    \\node [latent] (y1wt-1)  {$y_{1,t-1}^w$};
    \\node [latent] (y1wt)    [right =5em of y1wt-1]  {$y_{1,t}^w$};
    \\node [latent] (delta1t)   [right =2.5em of y1wt]  {$\\delta_{{y_1},t}$};

    % draw nodes
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
      VARmodel$q, "}
    \\node [latent] (y\\i wt-1) [below = 2.5em of y\\lasti wt-1] {$y_{\\i ,t-1}^w$};

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    VARmodel$q, "}
    \\node [latent] (y\\i wt) [below = 2.5em of y\\lasti wt] {$y_{\\i ,t}^w$};

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    VARmodel$q, "}
    \\node [latent] (delta\\i t) [below = 2.5em of delta\\lasti t] {$\\delta_{{y_\\i} ,t}$};

    % draw paths from residuals to yiwt
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
    VARmodel$q, "}
    \\draw [path] (delta\\i t) to node [] {} (y\\i wt);
    ")
    # WORK HERE FABI
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
              "\\draw [path", ifelse(
                # decorate with dot on path if parameter is random
                all_phis[all_phis$Param == paste0("phi_", i, j), "isRandom"] == 1,
                ", postaction = random]", "]"
              ),
              "  (y", i, "wt-1)  to node  []  {$\\phi_{y_", i, "}$}  (y", i, "wt);\n"
            )
          )
        } else {
          # for cross-lagged paths
          phi_vec <- c(
            phi_vec, paste0(
              "\\draw [path", ifelse(
                # decorate with dot on path if parameter is random
                all_phis[all_phis$Param == paste0("phi_", i, j), "isRandom"] == 1,
                ", postaction = random_cl]", "]"
              ),
              "  (y", i, "wt-1)  to node  [pos = .2]  {$\\phi_{y_{", i, j, "}}$}  (y", j, "wt);\n"
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
          "\\draw [var", ifelse(
            # decorate with dot on path if parameter is random
            all_sigmas[all_sigmas$Param == paste0("ln.sigma2_", i), "isRandom"] == 1,
            ", postaction = random]", "]"
          ),
          "  (delta", i, "t.120)   to node [] {$\\pi_{y_", i, "}$} (delta", i, "t.60);\n"
        )
      )
    }
    # do the same for innovation covariances
    # all_psis <- model[grepl("ln.sigma_", model$Param) & grepl("Fix", model$Type), ]
    # n_psi <- nrow(all_psis)
    psi_vec <- c()
    if (VARmodel$q == 1) {
      psi_vec <- paste0(
        "\\draw [cov", ifelse(
          # decorate with dot on path if parameter is random
          all_psis[all_psis$Param == "ln.sigma_12", "isRandom"] == 1,
          ", postaction = random]", "]"
        ),
        "  (delta1t.0)   to node [] {$\\pi_{y_{12}}$} (delta2t.0);\n"
      )
    } else {
      for (i in 1:(VARmodel$q - 1)) {
        for (j in (i + 1):VARmodel$q) {
          psi_vec <- c(
            psi_vec, paste0(
              "\\draw [cov", ifelse(
                # decorate with dot on path if parameter is random
                all_psis[all_psis$Param == paste0("ln.sigma_", i, j), "isRandom"] == 1,
                ", postaction = random]", "]"
              ),
              "  (delta", i, "t.0)   to node [] {$\\pi_{y_{", i, j, "}}$} (delta", j, "t.0);\n"
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
    begin_figure,
    wm_caption,
    begin_tikz,
    within_model,
    end_tikz,
    end_figure,
    sep = "\n"
  )

  # paste within_model to markdown
  cat(within_model, file = "pathmodel.rmd", append = TRUE)


  # between-model #############################################################

  # # between-model caption
  # bm_caption <- "\\caption*{Between-model.}"
  #
  # # draw nodes for q time-series constructs
  # if (VARmodel$q == 1) {
  #   bm <- "
  #   % draw between-level structural model
  #   \\node [latent] (mu_y1) {$\\mu_{y_1,t}$};
  #   \\node [latent] (phi_y1) [right = 1.5em of mu_y1] {$\\phi_{y_1}$};
  #   \\node [latent] (lnsigma2_y1) [right = 1.5em of phi_y1] {$log(\\pi_{y_1})$};
  #   "
  #   # draw (co-)variances conditional on existing correlations
  #   if (nrow(model[model$Param == "r_mu_1.phi_11"]) > 0) {
  #     r_mu_1.phi_11 <- "\\draw [cov] (mu_y1.south) to node [] {} (phi_y1.south);"
  #   }
  #   if (nrow(model[model$Param == "r_mu_1.ln.sigma2_1"]) > 0) {
  #     r_mu_1.ln.sigma2_1 <- "\\draw [cov] (mu_y1.south) to node [] {} (lnsigma2_y1.south);"
  #   }
  #   if (nrow(model[model$Param == "r_phi_11.ln.sigma2_1"]) > 0) {
  #     r_phi_11.ln.sigma2_1 <- "\\draw [cov] (phi_y1.south) to node [] {} (lnsigma2_y1.south);"
  #   }
  # }
  #
  # # paste together
  # between_model <- paste(bm, r_mu_1.phi_11, r_mu_1.ln.sigma2_1, r_phi_11.ln.sigma2_1, sep = "\n")
  #
  # # paste together
  # between_model <- paste(
  #   begin_figure,
  #   bm_caption,
  #   begin_tikz,
  #   between_model,
  #   end_tikz,
  #   end_figure,
  #   sep = "\n"
  # )
  #
  # # paste within_model to markdown
  # cat(between_model, file = "pathmodel.rmd", append = TRUE)




  # render markdown input #####################################################

  # render markdown
  rmarkdown::render(
    input = "pathmodel.rmd"
  )
}
