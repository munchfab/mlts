test_that("VARmodelBuild returns correct number of parameters", {
  # test number of returned model parameters for up to 10 constructs
  for (i in 1:10) {
    test_model <- suppressMessages(VARmodelBuild(q = i))
    n_within_pars <- 2*i + i^2 + (i * (i - 1)) / 2
    n_between_pars <- n_within_pars + (n_within_pars * (n_within_pars - 1)) / 2
    # number of within parameters
    expect_equal(
      nrow(test_model[
        test_model$Model == "Structural" & test_model$Level == "Within",
      ]),
      n_within_pars
    )
    # number of between parameters
    expect_equal(
      nrow(test_model[
        test_model$Model == "Structural" & test_model$Level == "Between",
      ]),
      n_between_pars
    )
  }
})

test_that("VARmodelBuild throws an error if p != length(q)", {
  expect_error(VARmodelBuild(q = 2, p = c(2, 2, 2)))
})
