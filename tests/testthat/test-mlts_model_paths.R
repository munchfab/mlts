test_that("Common pathmodels are rendered without errors", {
  # 1 construct, 1 indicator, no lagged effects
  model <- mlts_model(q = 1)
  expect_no_error(
    mlts_model_paths(model, file = "./tests/testthat/pathmodel.pdf")
  )
  # 2 construct, 1 indicator, no lagged effects
  model <- mlts_model(q = 2)
  expect_no_error(
    mlts_model_paths(model, file = "./tests/testthat/pathmodel.pdf")
  )
  # 1 construct, 3 indicator, no lagged effects
  model <- mlts_model(q = 1, p = 3)
  expect_no_error(
    mlts_model_paths(model, file = "./tests/testthat/pathmodel.pdf")
  )
  # 2 construct, 3 indicator, no lagged effects
  model <- mlts_model(q = 2, p = c(3, 3))
  expect_no_error(
    mlts_model_paths(model, file = "./tests/testthat/pathmodel.pdf")
  )
  # 2 constructs, 1 + 3 indicators, AR2-effects
  model <- mlts_model(q = 2, p = c(1, 3), max_lag = 2)
  expect_no_error(
    mlts_model_paths(model, file = "./tests/testthat/pathmodel.pdf")
  )
  # 1 constructs, 3 indicators, with btw_factor = F
  model <- mlts_model(q = 1, p = 3)
  model <- mlts_model_measurement(model, q = 1,  p = 3, btw_factor = F)
  expect_no_error(
    mlts_model_paths(model, file = "./tests/testthat/pathmodel.pdf")
  )
  # 2 constructs, 3 indicators each, with btw_factor = F
  model <- mlts_model(q = 2, p = c(3, 3))
  model <- mlts_model_measurement(model, q = 2,  p = c(3, 3), btw_factor = F)
  expect_no_error(
    mlts_model_paths(model, file = "./tests/testthat/pathmodel.pdf")
  )
  # remove files
  files <- list.files("./tests/testthat/")
  file.remove(files[grep("pathmodel", files)])
})
