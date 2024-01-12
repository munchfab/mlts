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
  	cov/.style = {<->, thick, > = {stealth[]}, bend right = 30},
  	% style to add a circle in the middle of a path
    random/.style = {postaction = {decorate, decoration = {markings, mark = at position .5 with {\\draw[fill = black] circle[radius = 2pt];}}}},
  	]"

  # end tikz picture
  end_tikz <- "\\end{tikzpicture}\n"


  # decomposition #############################################################

  dc_caption <- "\\caption*{Decomposition.}"

  if (VARmodel$q == 1) {
    dc <- "
    \\node [manifest] (y1t) {$y_{1,t}$};
    \\node [latent]   (y1wt)    [above = 2.5em of y1t]  {$y_{1,t}^w$};
    \\node [latent]   (mu_y1)   [below = 2.5em of y1t]  {$\\mu_{y_1,t}$};

    % draw paths
    \\draw [path] (y1wt) to node [] {} (y1t);
    \\draw [path] (mu_y1) to node [] {} (y1t);"
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
    if (model[model$Param == "phi_11" & model$Type == "Fix effect", "isRandom"] == 1) {
      phi_11 <- "\\draw [path, postaction = random]  (y1wt-1)  to node [] {$\\phi_{y_1}$}  (y1wt);"
    } else {
      phi_11 <- "\\draw [path]  (y1wt-1)  to node [] {$\\phi_{y_1}$}  (y1wt);"
    }
    if (model[model$Param == "ln.sigma2_1" & model$Type == "Fix effect", "isRandom"] == 1) {
      ln.sigma2_1 <- "% draw (co-)variances
      \\draw [var, postaction = random]	  (delta1t.30)   to node [] {$\\pi_{y_1}$} (delta1t.330);"
    } else {
      ln.sigma2_1 <- "% draw (co-)variances
      \\draw [var]	  (delta1t.30)   to node [] {$\\pi_{y_1}$} (delta1t.330);"
    }
    # paste together
    within_model <- paste(wm, phi_11, ln.sigma2_1, sep = "\n")
  }

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

  # between-model caption
  bm_caption <- "\\caption*{Between-model.}"

  # draw nodes for q time-series constructs
  if (VARmodel$q == 1) {
    bm <- "
    % draw between-level structural model
    \\node [latent] (mu_y1) {$\\mu_{y_1,t}$};
    \\node [latent] (phi_y1) [right = 1.5em of mu_y1] {$\\phi_{y_1}$};
    \\node [latent] (lnsigma2_y1) [right = 1.5em of phi_y1] {$log(\\pi_{y_1})$};
    "
    # draw (co-)variances conditional on existing correlations
    if (nrow(model[model$Param == "r_mu_1.phi_11"]) > 0) {
      r_mu_1.phi_11 <- "\\draw [cov] (mu_y1.south) to node [] {} (phi_y1.south);"
    }
    if (nrow(model[model$Param == "r_mu_1.ln.sigma2_1"]) > 0) {
      r_mu_1.ln.sigma2_1 <- "\\draw [cov] (mu_y1.south) to node [] {} (lnsigma2_y1.south);"
    }
    if (nrow(model[model$Param == "r_phi_11.ln.sigma2_1"]) > 0) {
      r_phi_11.ln.sigma2_1 <- "\\draw [cov] (phi_y1.south) to node [] {} (lnsigma2_y1.south);"
    }
  }

  # paste together
  between_model <- paste(bm, r_mu_1.phi_11, r_mu_1.ln.sigma2_1, r_phi_11.ln.sigma2_1, sep = "\n")

  # paste together
  between_model <- paste(
    begin_figure,
    bm_caption,
    begin_tikz,
    between_model,
    end_tikz,
    end_figure,
    sep = "\n"
  )

  # paste within_model to markdown
  cat(between_model, file = "pathmodel.rmd", append = TRUE)




  # render markdown input #####################################################

  # render markdown
  rmarkdown::render(
    input = "pathmodel.rmd"
  )
}
