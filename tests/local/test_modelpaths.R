devtools::load_all()

# with random effects
VARmodel <- VARmodelBuild(q = 2, p = c(2, 2))
# VARmodel
# VARmodel[VARmodel$Model == "Measurement",]
# VARmodelEval(VARmodel)
VARmodelPaths(VARmodel = VARmodel)


VARmodel <-  VARmodelBuild(q = 2, RE.pred = c("x", "z"), out.pred = "y")
VARmodelPaths(VARmodel = VARmodel)



