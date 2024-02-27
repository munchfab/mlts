devtools::load_all()

# with random effects
VARmodel <- VARmodelBuild(q = 2, p = c(2, 2))
# VARmodel
# VARmodel[VARmodel$Model == "Measurement",]
# VARmodelEval(VARmodel)
VARmodelPaths(VARmodel = VARmodel)



devtools::load_all()
VARmodel <-  VARmodelBuild(q = 2, RE.pred = c("x"))

# VARmodel <- VARmodelConstraints(VARmodel, FEis0 = "b_phi(1)_22.ON.x")
VARmodelPaths(VARmodel = VARmodel)



VARmodel <-  VARmodelBuild(q = 2, out.pred = c("y", "y2"))
VARmodelformula(VARmodel = VARmodel)

devtools::load_all()
VARmodel <-  VARmodelBuild(q = 2, p = c(2, 2))
VARmodelformula(VARmodel = VARmodel)


devtools::load_all()
VARmodel <-  VARmodelBuild(q = 2, maxLag = 2, p = c(3, 2))
# VARmodelPaths(VARmodel = VARmodel)
VARmodelformula(VARmodel = VARmodel)

# extract model infos
infos <- VARmodelEval(VARmodel)

# extract model data frame
model <- VARmodel












