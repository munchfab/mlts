test_that("returned data is a data frame", {
  expect_s3_class(dsem_data, data.frame)
})
