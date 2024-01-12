

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

VARmodel <- list()

VARmodel <- list(VARmodel = VARmodeldata, q = 1)

#######################################
devtools::load_all()

VARmodelPaths(VARmodel = VARmodel)



if (file.exists("pathmodel.rmd")) {
  file.remove("pathmodel.rmd")
}

rmarkdown::draft(file = "pathmodel.rmd",
                 template = "pathmodel",
                 package = "dsemr",
                 edit = FALSE)

pathmodel <- paste0(
"\\begin{tikzpicture}[
  auto, > = latex, align=center,
	latent/.style = {circle, draw, thick, inner sep = 2pt, minimum width = \\Radius},
	intercept/.style = {regular polygon,regular polygon sides = 3, draw, thick, inner sep = 0pt, minimum size = \\Radius},
	manifest/.style = {rectangle, draw, thick, inner sep = 0pt, minimum size = 5mm},
	mean/.style = {regular polygon, regular polygon sides = 3, draw, thick, inner sep = 0pt, minimum size = 8mm},
	paths/.style = {arrows = ->, thick, > = {stealth[]}},
	error/.style = {circle, draw = none, fill = none, thick, inner sep = 0pt, minimum size = 5mm},
	var/.style = {<->, thick, > = {stealth[]}, bend right = 270, looseness = 2},
	cov/.style = {<->, thick, > = {stealth[]}, bend right = 270},
	% style to add a circle in the middle of a path
  random/.style = {postaction = {decorate, decoration = {markings, mark = at position .5 with {\\draw[fill = black] circle[radius = 2pt];}}}},
	]

  % draw within-level structural model
  \\node [latent] (y1wt-1)                            {$y_{1,t-1}^w$};
  \\node [latent] (y1wt)    [right =5em of y1wt-1]    {$y_{1,t}^w$};
  \\node [latent] (eps1t)   [right =2.5em of y1wt]    {$\\delta_{{y_1},t}$};

  % draw paths
  \\draw [paths, postaction = random]  (y1wt-1)  to node [] {$\\phi_{y_1}$}  (y1wt);
  \\draw [paths]  (eps1t)   to node [] {}           (y1wt);

  % draw (co-)variances
  \\draw [var]	  (eps1t.30)   to node [] {} (eps1t.330);
  ",
  # if (VARmodel[VARmodel$Param == "mu_1" & VARmodel$Type == "Fix effect", "isRandom"] == 1) {
  #   "\\draw [circle] (0, 0) -- (2, 0);"
  # },

"\\end{tikzpicture}
"
)


cat(pathmodel, file = "model.rmd", append = TRUE)

rmarkdown::render(
  input = "model.rmd"
)




### Test model fitting function
#### 1. simulated data
load("data/var1_sample_data.Rdata", verbose = T)
VARresult = VARfit(VARmodel = VARmodel, data = dat$WITHIN, ts.ind = "X1", iter = 1500)

# get summary output
monitor(VARresult$stanfit, digits_summary = 3)
