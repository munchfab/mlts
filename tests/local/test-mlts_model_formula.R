test_that("Common model formulas are rendered without errors", {
  path <- "tests/local"
  file <- "tests/local/test_formula.pdf"
  # 1 construct, 1 indicator, no lagged effects
  model <- mlts_model(q = 1)
  expect_no_error(
    mlts_model_formula(model, file = file)
  )
  # 2 construct, 1 indicator, no lagged effects
  model <- mlts_model(q = 2)
  expect_no_error(
    mlts_model_formula(model, file = file)
  )
  # 1 construct, 3 indicator, no lagged effects
  model <- mlts_model(q = 1, p = 3)
  expect_no_error(
    mlts_model_formula(model, file = file)
  )
  # 2 construct, 3 indicator, no lagged effects
  model <- mlts_model(q = 2, p = c(3, 3))
  expect_no_error(
    mlts_model_formula(model, file = file)
  )
  # 2 constructs, 1 + 3 indicators, AR2-effects
  model <- mlts_model(q = 2, p = c(1, 3), max_lag = 2)
  expect_no_error(
    mlts_model_formula(model, file = file)
  )
  # 2 constructs, 2 + 3 indicators, AR3-effects
  model <- mlts_model(q = 2, p = c(2, 3), max_lag = 3)
  expect_no_error(
    mlts_model_formula(model, file = file)
  )
  # 1 constructs, 3 indicators, with btw_factor = F
  model <- mlts_model(q = 1, p = 3)
  model <- mlts_model_measurement(model, q = 1,  p = 3, btw_factor = F)
  expect_no_error(
    mlts_model_formula(model, file = file)
  )
  # 2 constructs, 3 indicators each, with btw_factor = F
  model <- mlts_model(q = 2, p = c(3, 3))
  model <- mlts_model_measurement(model, q = 2,  p = c(3, 3), btw_factor = F)
  expect_no_error(
    mlts_model_formula(model, file = file)
  )
  # 1 construct, 1 indicator, phi fixed
  model <- mlts_model(q = 1)
  model <- mlts_model_constraint(model, fix_dynamics = T)
  expect_no_error(
    mlts_model_formula(model, file = file)
  )
  # 2 construct, 3 indicator each, ln(sigma) fixed
  # doesnt work?
  model <- mlts_model(q = 2, p = c(3, 3))
  model <- mlts_model_constraint(model, fix_inno_vars = T)
  expect_no_error(
    mlts_model_formula(model, file = file)
  )
  # 3 construct, 3 indicator each, innovation covariance fixed
  model <- mlts_model(q = 3, p = c(3, 3, 3))
  model <- mlts_model_constraint(model, fix_inno_covs = T)
  expect_no_error(
    mlts_model_formula(model, file = file)
  )
  # 2 construct, 3 indicator each, AR-effect fixed
  model <- mlts_model(q = 2, p = c(3, 3))
  model <- mlts_model_constraint(model, fixef_zero = c("phi(1)_21"))
  expect_no_error(
    mlts_model_formula(model, file = file)
  )
  # 2 construct, 3 indicator each, some random effects fixed
  model <- mlts_model(q = 2, p = c(3, 3))
  model <- mlts_model_constraint(model, ranef_zero = c("phi(1)_21", "phi(1)_12"))
  expect_no_error(
    mlts_model_formula(model, file = file)
  )
  # 3 construct, 1 indicator each, all innovations independent?
  model <- mlts_model(q = 2)
  model <- mlts_model_constraint(model, inno_covs_zero = T)
  expect_no_error(
    mlts_model_formula(model, file = file)
  )
  # remove files
  files <- list.files(path, full.names = TRUE)
  file.remove(files[grep("test_formula", files)])
})

