devtools::load_all()

## AR(1) Models ================================================================
# 1. Random effect model ----------------------------------------------------
### Build general model
(VARmodeldata = VARmodelBuild(q = 1))

# create artificial data
ar1_data <- sim_data_AR1_man(
  N = 100,
  TP = 70,
  mu = 10, phi = .3, logv = 2,
  sigma_mu = 1.4, sigma_phi = .2, sigma_logv = .7,
  mu_ec = 5, sigma_ec = 2,
  cor_mu_phi = .3, cor_mu_logv = .3, cor_phi_logv = .3,
  cor_mu_ec = .3, cor_phi_ec = .3, cor_logv_ec = .3,
  seed = 1234
)

data <- ARprepare(VARmodel = VARmodel, data = ar1_data)

ar1_data_stan <- create_stan_data(
  data = ar1_data,
  y = "y", id = "id", beep = "TP",
  miss_handling = "remove"
)

(VARmodeldata = VARmodelBuild(q = 2))

VARmodel <- list()

VARmodel <- list(VARmodel = VARmodeldata, q = 2)




#######################################
devtools::load_all()

# with random effects
VARmodeldata = VARmodelBuild(q = 3)
VARmodel <- list()
VARmodel <- list(VARmodel = VARmodeldata, q = 3)

VARmodelPaths(VARmodel = VARmodel)




devtools::load_all()
# with random effects
VARmodeldata = VARmodelBuild(q = 2)
# VARmodeldata = VARmodelConstraints(VARmodel = VARmodeldata, FEis0 = c("phi_11"))
VARmodel <- list()
VARmodel <- list(VARmodel = VARmodeldata, q = 2)




VARmodelformula(VARmodel = VARmodel)




# render markdown
rmarkdown::render(
  input = "pathmodel.rmd"
)


# without random effects
VARmodeldata = VARmodelBuild(q = 1)
VARmodeldata[VARmodeldata$Param == "phi_11" & VARmodeldata$Type == "Fix effect", "isRandom"] <- 0
VARmodeldata[VARmodeldata$Param == "ln.sigma2_1" & VARmodeldata$Type == "Fix effect", "isRandom"] <- 0
VARmodel <- list()
VARmodel <- list(VARmodel = VARmodeldata, q = 1)

VARmodelPaths(VARmodel = VARmodel)



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

  dc <- paste0(
    "
    \\node [manifest] (y1t) {$y_{1,t}$};
    \\node [latent]   (y1wt)    [above = 2.5em of y1t]  {$y_{1,t}^w$};
    \\node [latent]   (mu_y1)   [below = 2.5em of y1t]  {$\\mu_{y_1,t}$};

    \\foreach \\i in {1, ...,",  VARmodel$q, "}
    \\node [manifest] (y\\it) [right = 2.5em of y\\it] {$y_{1,t}$};

    % draw paths
    \\draw [path] (y1wt) to node [] {} (y1t);
    \\draw [path] (mu_y1) to node [] {} (y1t);"
  )

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


# render markdown
rmarkdown::render(
  input = "pathmodel.rmd"
)



"\begin{tikzpicture}[
  auto, > = latex, align=center,
  latent/.style = {circle, draw, thick, inner sep = 2pt, minimum width = \Radius},
  manifest/.style = {rectangle, draw, thick, inner sep = 0pt, minimum size = \Radius/2},
  intercept/.style = {regular polygon,regular polygon sides = 3, draw, thick, inner sep = 0pt, minimum size = \Radius},
  mean/.style = {regular polygon, regular polygon sides = 3, draw, thick, inner sep = 0pt, minimum size = 8mm},
  path/.style = {arrows = ->, thick, > = {stealth[]}},
  error/.style = {circle, draw = none, fill = none, thick, inner sep = 0pt, minimum size = 5mm},
  var/.style = {<->, thick, > = {stealth[]}, bend right = 270, looseness = 2},
  cov/.style = {<->, thick, > = {stealth[]}, bend right = 30},
  % style to add a circle in the middle of a path
  random/.style = {postaction = {decorate, decoration = {markings, mark = at position .5 with {\draw[fill = black] circle[radius = 2pt];}}}},
]

%\node [manifest] (y1t) {$y_{1,t}$};

%works?
  %\foreach \x/\xtext in {0,...,3}
%\node [manifest] (\xtext) at (\x,0) {$\xtext$};


% works for arrows
%\foreach \x [remember=\x as \lastx (initially 1)] in {1,...,3} {\lastx$\to$\x, };

% initiate first
\node [manifest] (0) {0};
\foreach \x [remember=\x as \lastx (initially 0)] in {1,...,3} \node [manifest] (\x) [right =2.5em of \lastx] {$\x$};



%foreach \i/\itext in {0,...,3}
%\node [manifest] (\itext) [right = 2.5em of \itext-1] {$\itext$};

\end{tikzpicture}

\end{figure}"
