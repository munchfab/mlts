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
               M           = 0.9,
               Mdn         = 0.75,
               SD          = 1.062,
               Skew        = 0.272,
               Kurtosis    = 0.605,
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
                 M           = c(0.647, -0.026),
                 Mdn         = c(0.601, -0.127),
                 SD          = c(1.254, 1.286),
                 Skew        = c(0.232, 0.183),
                 Kurtosis    = c(-0.345, -0.045),
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
                 M           = c(0.338, 0.764, 2.229),
                 Mdn         = c(0.409, 0.720, 2.217),
                 SD          = c(1.161, 0.999, 1.025),
                 Skew        = c(0.112, 0.475, 0.047),
                 Kurtosis    = c(-0.111, -0.250, -0.452),
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
                 M           = c(0.626, 0.981, 2.514, 0.641, 1.934, 1.105),
                 Mdn         = c(0.587, 0.954, 2.441, 0.554, 1.864, 1.084),
                 SD          = c(1.132, 0.829, 0.967, 1.261, 0.985, 1.127),
                 Skew        = c(0.232, -0.028,  0.177,  0.098,  0.187,  0.038),
                 Kurtosis    = c(-0.378,  0.116, -0.423,  0.153, -0.275,  0.259),
                 floor_pc    = c(1, 1, 1, 1, 1, 1),
                 ceiling_pc  = c(1, 1, 1, 1, 1, 1)
               ))
  txt <- "
group ts_var  N N_comp miss_pc N_lag1_comp     M   Mdn    SD   Skew Kurtosis floor_pc ceiling_pc
1      1     Y1 10     10       0           9 0.676 0.926 1.558 -0.515   -1.334       10          0
2      2     Y1 10     10       0           9 0.908 0.793 0.748  0.246   -1.466        0          0
3      3     Y1 10     10       0           9 0.337 0.376 0.471 -0.258   -0.637        0          0
4      4     Y1 10     10       0           9 0.137 0.333 0.633 -0.570   -0.404        0          0
5      5     Y1 10     10       0           9 0.128 0.286 0.549 -1.674    3.459        0          0
6      6     Y1 10     10       0           9 0.998 1.018 0.599  0.431    0.382        0          0
7      7     Y1 10     10       0           9 0.950 1.017 0.833 -1.380    2.572        0          0
8      8     Y1 10     10       0           9 1.701 1.539 1.429  0.665    0.201        0         10
9      9     Y1 10     10       0           9 1.394 1.520 1.109 -0.481   -0.928        0          0
10    10     Y1 10     10       0           9 1.770 2.148 0.783 -0.727   -0.963        0          0
"

  df1_group <- read.table(text = txt, header = TRUE, row.names = 1)
  row.names(df1_group) <- NULL
  expect_equal(describe_ts(data = initial_simData1$data, group = "ID", ts = "Y1"), df1_group)
})
