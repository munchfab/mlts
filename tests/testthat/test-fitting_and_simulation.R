test_that("fitting", {
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
                        data = initial_simData1$data,
                        id = "ID",
                        ts = "Y1",
                        chains = 2,
                        iter = 5,
                        cores = 1,
                        seed = 5)$pop.pars.summary

  new_fit2 <- mlts_fit(model = m2,
                        data = initial_simData2$data,
                        id = "ID",
                        ts = c("Y1","Y2"),
                        chains = 2,
                        seed = 5,
                        iter = 5,
                        cores = 1)$pop.pars.summary

  new_fit3 <- mlts_fit(model = m3,
                       data = initial_simData3$data,
                       id = "ID",
                       ts = c("Y1.1","Y1.2", "Y1.3"),
                       chains = 2,
                       seed = 5,
                       iter = 5,
                       cores = 1)$pop.pars.summary

  new_fit4 <- mlts_fit(model = m4,
                       data = initial_simData4$data,
                       id = "ID",
                       ts = c("Y1.1", "Y1.2", "Y1.3",
                              "Y2.1", "Y2.2", "Y2.3"),
                       chains = 2,
                       seed = 5,
                       iter = 5,
                       cores = 1)$pop.pars.summary


  expect_equal(new_fit1[, "mean"], c( -0.045,  0.423,  0.870,  2.994,  0.292,
                                      4.488, 0.115, -0.704, -0.209),
               tolerance = 1e-5)
  expect_equal(new_fit2[, "mean"], c(-0.099, -1.465, -1.605, -0.532,  0.809,
                                     -1.378, 0.184,  1.793,  0.740,  1.394,
                                     0.544,  0.063),
               tolerance = 1e-5)
  expect_equal(new_fit3[, "mean"], c( -0.125, -0.582,  0.995,  1.682,  0.000,
                                      0.348, -1.194,  1.000,  0.094, -0.198,
                                      1.000, -0.054, -1.183,  0.000,  1.467,
                                      0.609,  1.744, 1.279, 2.808),
               tolerance = 1e-5)
  expect_equal(new_fit4[, "mean"], c(-0.080, -1.461, -1.667, -0.547,  0.844,
                                     -1.414,  0.172,  1.730,  0.578,  1.392,
                                      0.542,  0.057,  0.000, -1.539, -0.177,
                                      0.000, -0.367,  0.440,  1.000,  1.339,
                                      0.288,  1.000, -0.886,  0.206,  1.000,
                                     -0.036,  0.531,  1.000,  0.869,  0.510,
                                      0.000,  2.875,  3.778,  0.000,  2.414,
                                      2.752,  0.774,  2.546,  2.898,  2.038,
                                      0.991,  2.636),
               tolerance = 1e-5)

})


# testing
test_that("simulation", {
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
               c(0.80,  0.15, -0.30,  1.10,  0.15,  0.25, -0.20,  0.05, -0.20))
  expect_equal(new_simData2$model[, "true.val"],
               c(0.80,  0.30,  0.15, -0.20,  0.05,  0.15,  0.75,  0.75, -0.15,
                 0.70,  1.20, -0.20))
  expect_equal(new_simData3$model[, "true.val"],
               c(0.40, 0.20, 0.75, 0.90, 0.00, 0.50, 2.00, 1.00, 0.90, 0.80,
                 0.00, 0.15, 0.15, 1.00, 0.70, 0.75, 0.20, 0.20, 0.25))
  expect_equal(new_simData4$model[, "true.val"],
               c(0.80,  0.80,  0.25, -0.05, -0.05,  0.25,  0.75,  0.75, -0.15,
                 0.70,  0.70,  0.15,  0.00,  0.50,  2.00,  0.00,  1.50,  0.50,
                 1.00,  0.80,  0.80,  1.00,  0.70,  0.90,  0.00,  0.15,  0.15,
                 0.00,  0.15,  0.15,  1.00,  0.75,  0.90,  1.00,  0.80,  0.75,
                 0.15,  0.20,  0.20,  0.20,  0.20,  0.25))
})


