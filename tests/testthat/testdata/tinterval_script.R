library(mlts)

# load data
dat = read.csv(file = "data/Wenzel_Brose_2022_dataset3.csv")
# sub = dat[dat$id == min(dat$id),c("id", "time", "dayno", "beep", "pos")]
# save(sub, file = "tinterval_data.rda")
load(file = "./tests/testthat/testdata/tinterval_data.rda")


# plot data of first individual
sub = sub[sub$time < 10,]          # first 10 hours
plot(y = sub$pos, x = sub$time)
lines(y = sub$pos, x = sub$time)
sub

# illustration of time-grid
sub.nas = prepare_data(sub, id = "id", ts = "pos", time = "time", tinterval = 1)
sub.nas[,c("id","time", "int_time", "pos")]
sub.nas


# illustration of time-grid
sub.nas2 = prepare_data(sub, id = "id", ts = "pos", time = "time", tinterval = 1)
sub.nas2[,c("id","time", "int_time", "pos")]
sub.nas2


# create some data for example
data <- data.frame(
  id = rep(c(1, 2), each = 4),
  time = c(0, 3, 4, 6,
           1, 4, 5, 7)
)

# create missings to approximate continuous time
create_missings(
  data = data, id = "id", time = "time",
  tinterval = 1 # use time interval of 1 minute
)
