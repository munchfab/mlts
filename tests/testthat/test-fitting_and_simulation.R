test_that("Das Modell schätzt die gleichen Werte wie immer ", {
  skip_on_cran()

  m1 <- mlts_model(q = 1)
  m2 <- mlts_model(q = 2, fix_dynamics = TRUE, fix_inno_vars = TRUE)
  m3 <- mlts_model(q = 1, p = 3, fix_dynamics = TRUE, fix_inno_vars = TRUE)
  m4 <- mlts_model(q = 2, p = c(3, 3), fix_dynamics = TRUE, fix_inno_vars = TRUE)
  # Lade die Referenzdaten
  initial_simData1  <- readRDS(testthat::test_path("testdata",
                                                   "simData1.rds"))
  initial_simData2  <- readRDS(testthat::test_path("testdata",
                                                   "simData2.rds"))
  initial_simData3  <- readRDS(testthat::test_path("testdata",
                                                   "simData3.rds"))
  initial_simData4  <- readRDS(testthat::test_path("testdata",
                                                   "simData4.rds"))

  # neues fitting
  new_fit1 <- mlts_fit(model = m1,
                        data = initial_simData1,
                        id = "ID",
                        ts = "Y1",
                        chains = 1,
                        iter = 5,
                        cores = 1,
                        seed = 5)$pop.pars.summary

  new_fit2 <- mlts_fit(model = m2,
                        data = initial_simData2,
                        id = "ID",
                        ts = c("Y1","Y2"),
                        chains = 1,
                        seed = 5,
                        iter = 5,
                        cores = 1)$pop.pars.summary

  new_fit3 <- mlts_fit(model = m3,
                       data = initial_simData3,
                       id = "ID",
                       ts = c("Y1.1","Y1.2", "Y1.3"),
                       chains = 1,
                       seed = 5,
                       iter = 5,
                       cores = 1)$pop.pars.summary

  new_fit4 <- mlts_fit(model = m4,
                       data = initial_simData4,
                       id = "ID",
                       ts = c("Y1.1", "Y1.2", "Y1.3",
                              "Y2.1", "Y2.2", "Y2.3"),
                       chains = 1,
                       seed = 5,
                       iter = 5,
                       cores = 1)$pop.pars.summary


  expect_equal(new_fit1[, "mean"], c(-0.025,  0.072,  0.662,  1.212,
                                      0.215,  7.029, -0.632, -0.633, -0.026),
               tolerance = 1e-5)
  expect_equal(new_fit2[, "mean"], c(-1.738, -1.366, -1.838,  0.025,  0.176,
                                     -1.502, -0.408,  3.181,  1.100,  0.514,
                                      0.915, 0.896),
               tolerance = 1e-5)
  expect_equal(new_fit3[, "mean"], c( 0.373, -0.208,  1.567,  0.522,  0.000,
                                     -0.464, -1.257,  1.000,  1.112, -1.071,
                                      1.000, -0.035, -1.085,  0.000,  0.683,
                                      1.020,  2.850,  1.969,  4.661),
               tolerance = 1e-5)
  expect_equal(new_fit4[, "mean"], c(-1.700, -1.359, -1.847,  0.019,  0.191,
                                     -1.585, -0.554,  3.193,  0.800,  0.514,
                                      0.910,  0.883,  0.000, -1.255, -1.768,
                                      0.000, -0.826,  1.474,  1.000,  1.999,
                                     -0.903,  1.000, -0.756, -1.068,  1.000,
                                      0.004,  0.199,  1.000,  0.668, -0.329,
                                      0.000,  0.691,  1.579,  0.000,  0.155,
                                      1.539,  1.159,  4.858,  1.227,  3.087,
                                      1.423,  5.064),
               tolerance = 1e-5)

})


# testing
# vermutlich nicht sehr sinnvoll, wenn default = T, weil dann sowieso immer
# gleiche pop pars ausgewählt werden - habe es trotzdem mal hinzugefügt
test_that("Neue simulierte Daten mit gleichem Seed die gleichen pop.pars haben", {
  skip_on_cran()

  m1 <- mlts_model(q = 1)
  m2 <- mlts_model(q = 2, fix_dynamics = TRUE, fix_inno_vars = TRUE)
  m3 <- mlts_model(q = 1, p = 3, fix_dynamics = TRUE, fix_inno_vars = TRUE)
  m4 <- mlts_model(q = 2, p = c(3, 3), fix_dynamics = TRUE, fix_inno_vars = TRUE)
  # neue Simulation
  new_simData1 = mlts_sim(model = m1,
                          N = 10,
                          TP = 10,
                          burn.in = 1,
                          default = TRUE,
                          seed = 5)
  new_simData2 = mlts_sim(model = m2,
                          N = 10,
                          TP = 10,
                          burn.in = 1,
                          default = TRUE,
                          seed = 5)
  new_simData3 = mlts_sim(model = m3,
                          N = 10,
                          TP = 10,
                          burn.in = 1,
                          default = TRUE,
                          seed = 5)
  new_simData4 = mlts_sim(model = m4,
                          N = 10,
                          TP = 10,
                          burn.in = 1,
                          default = TRUE,
                          seed = 5)

  expect_equal(new_simData1$model[, "true.val"],
               c(0.800,  0.300, -0.300,  0.700,  0.150,  0.250, -0.175,  0.050,
                 0.125))
  expect_equal(new_simData2$model[, "true.val"],
               c(0.80,  0.30,  0.25, -0.15,  0.00,  0.15,  0.75,  0.75, -0.15,
                 0.90,  1.20, -0.20))
  expect_equal(new_simData3$model[, "true.val"],
               c(0.40, 0.15, 0.75, 0.80, 0.00, 0.50, 2.00, 1.00, 0.90, 0.80,
                 0.00 ,0.15, 0.15, 1.00, 0.70, 0.75, 0.20, 0.20, 0.25))
  expect_equal(new_simData4$model[, "true.val"],
               c( 0.800,  1.000,  0.200,  0.000, -0.150,  0.150,  0.750,  0.750,
                 -0.150,  0.800,  1.200,  0.075,  0.000,  0.500,  2.000,  0.000,
                  1.500,  0.500,  1.000,  0.800,  0.800,  1.000,  0.700,  0.900,
                  0.000,  0.150,  0.150,  0.000,  0.150,  0.150,  1.000,  0.750,
                  0.900,  1.000,  0.800,  0.750,  0.150,  0.200,  0.200,  0.200,
                  0.200,  0.250))
})


