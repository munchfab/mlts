test_that("Missings are inserted correctly (at least same as in MPlus)", {
  expect_equal(2 * 2, 4)
})

mplus_data <- read.table("./data/missing_imputation_30minutes1.dat",
                         header = F)

colnames(mplus_data) <- c("V1", "V2", "timecont_day", "timeint", "V5",
                          "id", "int_time")

mplus_data <- mplus_data[order(mplus_data$id), ]
mplus_data <- as.data.frame(apply(
  mplus_data, 2,FUN = function(x) ifelse(x == "*", NA, as.numeric(x))
))

mplus15 <- mplus_data[, c("id", "timecont_day", "timeint", "int_time", "V1")]
mplus15$num_id <- as.numeric(as.factor(mplus15$id))

summary(mplus15)

library(dplyr)
library(tidyr)

mplus <- mplus15 %>%
  group_by(id) %>%
  fill(timecont_day, .direction = "down") %>%
  mutate(
    diff = timecont_day - lag(timecont_day),
    timeint = ifelse(diff < 0, timecont_day, timeint),
    timecont = ifelse(is.na(diff),
                      timecont_day, lag(timecont_day) + timeint),
  ) %>%
  # fill(timecont, .direction = "down") %>%
  mutate(
    timecont2 = ifelse(diff < 0,
                       cumsum(coalesce(timecont, 0)) + timecont*0, timecont)
  )


mplus <- mplus15 %>%
  group_by(id) %>%
  fill(timecont_day, .direction = "down") %>%
  mutate(
    diff = timecont_day - lag(timecont_day),
    timeint = ifelse(int_time == 1 | diff < 0, timecont_day, timeint)
  ) %>%
  mutate(
    timecont = cumsum(coalesce(timeint, 0)) + timeint*0, timeint
  )



load("./data/example_data.Rdata")
load("./data/mplus15.rda")


d <- create_missings(dat, tinterval = 1000, id = "UserID", time = "timecont")


devtools::load_all()
d2 <- create_missings2(data = dat, tinterval = 30, id = "UserID", time = "timecont")
warnings()

all(d$timecont == mplus15$timecont, na.rm = T)

all(d2$int_time == mplus$int_time, na.rm = T)

which(d$timecont != mplus$timecont)
which(d2$timecont != mplus$timecont)

d[105, ]
mplus[105, ]
d2[99, ]
mplus[99, ]
