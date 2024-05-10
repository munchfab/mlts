test_that("Missings are inserted correctly (at least same as in MPlus)", {
  # load example data
  mplus_raw <- testthat::test_path("testdata", "mplus_raw.rda")
  load(mplus_raw)
  # mplus dataset with 30 minutes tinterval
  mplus30 <- testthat::test_path("testdata", "mplus30.rda")
  load(mplus30)
  # load("mplus30.rda")
  data30 <- create_missings(
    data = mplus_raw, tinterval = 30, id = "UserID", time = "timecont"
  )
  expect_equal(
    object = is.na(data30$var1),
    expected = is.na(mplus30$var1)
  )
  expect_equal(
    object = data30$int_time,
    expected = mplus30$int_time
  )
  # mplus dataset with 15 minutes tinterval
  mplus15 <- testthat::test_path("testdata", "mplus15.rda")
  load(mplus15)
  data15 <- create_missings(
    data = mplus_raw, tinterval = 15, id = "UserID", time = "timecont"
  )
  expect_equal(
    object = is.na(data15$var1),
    expected = is.na(mplus15$var1)
  )
  expect_equal(
    object = data15$int_time,
    expected = mplus15$int_time
  )
  ##
  # slightly more complicated example
  mplus_raw2 <- testthat::test_path("testdata", "mplus_raw2.rda")
  load(mplus_raw2)
  mplus2 <- testthat::test_path("testdata", "mplus2.rda")
  load(mplus2)
  data2 <- create_missings(
    data = mplus_raw2, tinterval = 2, id = "id", time = "time_cont"
  )
  expect_equal(
    object = is.na(data2$var1),
    expected = is.na(mplus2$var1)
  )
  expect_equal(
    object = data2$int_time,
    expected = mplus2$int_time
  )
})

