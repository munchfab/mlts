test_that("calculating descriptives_ts works", {
  skip_on_cran()

  initial_simData1  <- readRDS(testthat::test_path("testdata",
                                                   "simData1.rds"))
  initial_simData2  <- readRDS(testthat::test_path("testdata",
                                                   "simData2.rds"))
  initial_simData3  <- readRDS(testthat::test_path("testdata",
                                                   "simData3.rds"))
  initial_simData4  <- readRDS(testthat::test_path("testdata",
                                                   "simData4.rds"))
  expect_equal(describe_ts(data = initial_simData1$data, ts = "Y1"),
               data.frame(
               ts_var      = "Y1",
               N           = 100,
               N_comp      = 100,
               miss_pc     = 0,
               N_lag1_comp = 99,
               M           = 0.946,
               Mdn         = 0.879,
               SD          = 1.209,
               Skew        = 0.283,
               Kurtosis    = 0.253,
               floor_pc    = 1,
               ceiling_pc  = 1
  ))
  expect_equal(describe_ts(data = initial_simData2$data, ts = c("Y1", "Y2")),
               data.frame(
                 ts_var      = c("Y1", "Y2"),
                 N           = c(100, 100),
                 N_comp      = c(100, 100),
                 miss_pc     = c(0, 0),
                 N_lag1_comp = c(99, 99),
                 M           = c(0.678, -0.026),
                 Mdn         = c(0.728, -0.138),
                 SD          = c(1.110, 1.279),
                 Skew        = c(0.120, 0.187),
                 Kurtosis    = c(-0.433, 0.041),
                 floor_pc    = c(1, 1),
                 ceiling_pc  = c(1, 1)
               ))
  expect_equal(describe_ts(data = initial_simData3$data, ts = c("Y1.1", "Y1.2", "Y1.3")),
               data.frame(
                 ts_var      = c("Y1.1", "Y1.2", "Y1.3"),
                 N           = c(100, 100, 100),
                 N_comp      = c(100, 100, 100),
                 miss_pc     = c(0, 0, 0),
                 N_lag1_comp = c(99, 99, 99),
                 M           = c(0.332, 0.759, 2.224),
                 Mdn         = c(0.323, 0.698, 2.189),
                 SD          = c(1.239, 1.073, 1.090),
                 Skew        = c(0.155, 0.503, 0.092),
                 Kurtosis    = c(-0.163, -0.280, -0.486),
                 floor_pc    = c(1, 1, 1),
                 ceiling_pc  = c(1, 1, 1)
               ))
  expect_equal(describe_ts(data = initial_simData4$data,
                           ts = c("Y1.1", "Y1.2", "Y1.3", "Y2.1", "Y2.2", "Y2.3")),
               data.frame(
                 ts_var      = c("Y1.1", "Y1.2", "Y1.3", "Y2.1", "Y2.2", "Y2.3"),
                 N           = c(100, 100, 100, 100, 100, 100),
                 N_comp      = c(100, 100, 100, 100, 100, 100),
                 miss_pc     = c(0, 0, 0, 0, 0, 0),
                 N_lag1_comp = c(99, 99, 99, 99, 99, 99),
                 M           = c(0.647, 0.997, 2.530, 0.613, 1.916, 1.078),
                 Mdn         = c(0.639, 0.988, 2.484, 0.506, 1.837, 1.077),
                 SD          = c(1.086, 0.790, 0.938, 0.964, 0.787, 0.821),
                 Skew        = c(0.132, -0.149,  0.079,  0.250,  0.149,  0.146),
                 Kurtosis    = c(-0.427,  0.204, -0.507,  0.235, -0.060,  0.066),
                 floor_pc    = c(1, 1, 1, 1, 1, 1),
                 ceiling_pc  = c(1, 1, 1, 1, 1, 1)
               ))
  txt <- "
group ts_var  N N_comp miss_pc N_lag1_comp      M    Mdn    SD   Skew Kurtosis floor_pc ceiling_pc
1      1     Y1 10     10       0           9  0.224  0.665 1.412 -0.577   -1.280       10          0
2      2     Y1 10     10       0           9  0.935  0.767 0.728  0.490   -1.058        0          0
3      3     Y1 10     10       0           9  0.154  0.203 0.522  0.068   -0.746        0          0
4      4     Y1 10     10       0           9  0.330  0.475 0.586 -0.858    0.631        0          0
5      5     Y1 10     10       0           9 -0.313 -0.196 0.567 -1.568    3.011        0          0
6      6     Y1 10     10       0           9  0.928  0.955 0.623  0.493    0.509        0          0
7      7     Y1 10     10       0           9  1.063  1.026 0.839 -1.234    2.409        0          0
8      8     Y1 10     10       0           9  2.032  2.030 1.327  0.636    0.052        0         10
9      9     Y1 10     10       0           9  1.700  1.799 1.152 -0.344   -1.586        0          0
10    10     Y1 10     10       0           9  2.402  2.700 0.810 -0.716   -0.529        0          0
"

  df1_group <- read.table(text = txt, header = TRUE, row.names = 1)
  row.names(df1_group) <- NULL
  expect_equal(describe_ts(data = initial_simData1$data, group = "ID", ts = "Y1"), df1_group)
})


