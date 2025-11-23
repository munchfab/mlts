test_that("Das Modell schätzt die gleichen Werte wie immer ", {
  skip_on_cran()
  m1 <- mlts_model(q = 1)

  # Lade die Referenzdaten
  path5  <- testthat::test_path("testdata", "simData0_t0_seed5.rds")
  path27 <- testthat::test_path("testdata", "simData0_t0_seed27.rds")
  initial_simData_seed5  <- readRDS(path5)
  initial_simData_seed27 <- readRDS(path27)

  # neues fitting
  new_fit_seed5 <- mlts_fit(model = m1, data = initial_simData_seed5,
                            id = "ID", ts = "Y1",
                            chains = 1, iter = 5,
                            cores = 1, seed = 5)$pop.pars.summary

  new_fit_seed27 <- mlts_fit(model = m1, data = initial_simData_seed27,
                             id = "ID", ts = "Y1",
                             chains = 1, iter = 5,
                             cores = 1, seed = 27)$pop.pars.summary

  expect_equal(new_fit_seed5[, "mean"], c(-0.025,  0.072,  0.662,  1.212,
                                          0.215,  7.029, -0.632, -0.633, -0.026),
               tolerance = 1e-5)
  expect_equal(new_fit_seed27[, "mean"], c(-1.605,  1.587,  1.311,  6.459,
                                           1.603,  6.767, -0.740,  0.626,  0.004),
               tolerance = 1e-5)
})


# testing
# vermutlich nicht sehr sinnvoll, wenn default = T, weil dann sowieso immer
# gleiche pop pars ausgewählt werden - habe es trotzdem mal hinzugefügt
test_that("Neue simulierte Daten mit gleichem Seed die gleichen pop.pars haben", {

  m1 <- mlts_model(q = 1)
  skip_on_cran()
  # neue Simulation
  new_simData_seed5 = mlts_sim(model = m1,
                               N = 10,
                               TP = 10,
                               burn.in = 1,
                               default = T,
                               seed = 5)

  new_simData_seed27 = mlts_sim(model = m1,
                                N = 10,
                                TP = 10,
                                burn.in = 1,
                                default = T,
                                seed = 27)

  expect_equal(new_simData_seed5$model[, "true.val"],
               c(0.800,  0.300, -0.300,  0.700,  0.150,  0.250, -0.175,  0.050,  0.125))
  expect_equal(new_simData_seed27$model[, "true.val"],
               c(0.800,  0.300, -0.300,  0.700,  0.150,  0.250, -0.175,  0.050,  0.125))
})
