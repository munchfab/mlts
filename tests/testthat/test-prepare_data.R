test_that("prepare_data() throws error if variables are missing", {
  # load example data
  mplus_raw <- testthat::test_path("testdata", "mplus_raw.rda")
  load(mplus_raw)
  # missing ts argument
  expect_error(
    prepare_data(mplus_raw, id = "UserID", tinterval = 1)
  )
  # missing time argument
  expect_error(
    prepare_data(mplus_raw, id = "UserID", tinterval = 1, ts = "var1")
  )
  # all set
  expect_no_error(
    prepare_data(mplus_raw, id = "UserID", tinterval = 1, ts = "var1",
                 time = "timecont")
  )
})
