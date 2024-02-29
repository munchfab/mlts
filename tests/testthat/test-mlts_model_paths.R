test_that("Common pathmodels are rendered without errors", {
  # 1 construct, 1 indicator, no lagged effects
  model <- mlts_model(q = 1)
  expect_no_error(
    mlts_model_paths(model)
  )
  # 2 construct, 3 indicator, no lagged effects
  model <- mlts_model(q = 1, p = 3)
  expect_no_error(
    mlts_model_paths(model)
  )
  # 2 constructs, 1 + 3 indicators, AR2-effects
  model <- mlts_model(q = 2, p = c(1, 3), max_lag = 2)
  expect_no_error(
    mlts_model_paths(model)
  )
})
