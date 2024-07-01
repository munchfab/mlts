#' Create Missings for Approximation of Continuous Time Dynamic Models
#' (new version)
#'
#' @param data An object of class `data.frame` (or one that can be coerced
#' to that class) containing data of all variables used in the model.
#' @param tinterval The step interval for approximation for a continuous time DSEM.
#' The smaller the step interval, the better the approximation.
#' @param id The variable in `data` that identifies the person or observational
#' unit (as character).
#' @param time The variable in `data` that contains the (continuous) time (as
#' string).
#' @param btw_vars The names of between-level variables in the data to be
#' added in newly created rows with NAs.
#'
#' @return A `data.frame` with missings imputed for use in \code{\link{mlts_fit}}.
#' @export
#'
#' @examples
#' # create some data for example
#' data <- data.frame(
#'   id = rep(c(1, 2), each = 4),
#'   time = c(0, 3, 4, 6,
#'            1, 4, 5, 7)
#' )
#'
#' # create missings to approximate continuous time process
#' create_missings(
#'   data = data, id = "id", time = "time",
#'   tinterval = 1 # use time interval of 1 minute
#' )
create_missings <- function(data, tinterval, id, time,
                            btw_vars = NULL) {

  # create empty list for imputed data frames
  imp_list <- list()
  # coerce data to data frame if necessary
  data <- as.data.frame(data)
  # create numeric id variable
  data$num_id <- as.numeric(as.factor(data[, id]))

  # remove between-level variables to add them later again
  btw_vars <- c(id, btw_vars)
  btw_data <- unique(data[, c("num_id", btw_vars)])
  data[, btw_vars] <- NULL

  for (i in 1:length(unique(data$num_id))) {

    # id
    num_id <- i

    # store continuous time variable for id == i and sort
    time_cont <- sort(data[data$num_id == i, time])

    # rescale time_cont
    # time_cont <- ceiling(time_cont / tinterval)
    # stop if time_cont is not strictly increasing
    if (any(diff(time_cont) <= 0)) {
      err_id <- btw_data[btw_data$num_id == i, id]
      # err_pos <- rownames(data[data$num_id == i, ])
      stop(
        time, " is not strictly increasing for id ", err_id,
        ". Please check entry ", which(diff(time_cont) == 0),
        " for id ", err_id, "."
      )
    }

    # seq in steps of time grid from min to max of time
    time_seq <- seq(min(time_cont), max(time_cont), by = tinterval)

    # indicate shift interval
    intervals <- findInterval(x = time_cont, vec = time_seq)
    max_shifts <- max(table(intervals))

    # print warning if max_shifts increases 4
    if (max_shifts > 4) {
      warn_id <- btw_data[btw_data$num_id == i, id]
      warning("Observations shifted ", max_shifts, " intervals for id ", warn_id)
    }

    # create empty vector of imputed time
    time_imp <- NA

    # check intervals for overlaps
    diff_intervals <- c(NA, diff(intervals))
    # store in temporary data frame
    temp <- data.frame(time_cont, intervals, diff_intervals)
    # search for overlapping intervals
    overlap <- which(temp$diff_intervals == 0)
    # search for nonoverlapping space between intervals
    nonoverlap <- which(temp$diff_intervals > 1)
    # create indicator variables for (non-)overlapping time intervals in temp
    temp[overlap, "overlap"] <- 1
    temp[nonoverlap, "nonoverlap"] <- 1
    # create new intervals_shifted variable from intervals
    temp$intervals_shifted <- temp$intervals
    # count number of times where observations are in same interval
    n_overlap <- length(overlap)
    # start while loop: while there is overlap, shift observations
    # to nearest empty interval
    while (n_overlap > 0) {
      # determine closest empty intervals
      empty_interval <- nonoverlap - overlap[1]
      empty_pos <- which.min(abs(nonoverlap - overlap[1]))
      # determine if shifts need to be up or down in data frame
      direction <- ifelse(
        empty_interval[empty_pos] > 0, "up", "down"
      )
      # shift according to direction
      if (direction == "down") {
        temp[nonoverlap[empty_pos]:(overlap[1] - 1), "intervals_shifted"] <-
          temp[nonoverlap[empty_pos]:(overlap[1] - 1), "intervals_shifted"] - 1
      } else {
        temp[overlap[1]:(nonoverlap[empty_pos] - 1), "intervals_shifted"] <-
          temp[overlap[1]:(nonoverlap[empty_pos] - 1), "intervals_shifted"] + 1
      }
      # calculate remaining overlap and start new if n_overlap > 0
      overlap <- which(diff(temp$intervals_shifted) == 0) + 1
      nonoverlap <- which(diff(temp$intervals_shifted) > 1) + 1
      nonoverlap[length(nonoverlap) + 1] <- dim(temp)[1] + 1
      n_overlap <- length(overlap)
    }

    for (k in 1:length(time_cont)) {
      time_imp[
        # determine in which interval between two numbers of time_seq
        # lies value [j] of time_cont
        temp$intervals_shifted[k]
      ] <- time_cont[k] # and fill that position with with value [j]
    }

    # create integer time variable for use in dsem analysis
    int_time <- seq_along(time_imp)
    # store data frame for id == i in imp_list[[i]]
    imp_list[[i]] <- data.frame(
      time_imp, int_time,
      num_id = rep(num_id, times = length(int_time))
    )
    colnames(imp_list[[i]]) <- c(time, "int_time", "num_id")

  }
  # reduce to data frame
  imp_data <- do.call("rbind", imp_list)

  # join with existing data set
  # with base R, somewhat slower than with dplyr
  new_data <- merge(x = imp_data, y = data, by = c("num_id", time), all = TRUE)
  # sort by num_id and int_time
  new_data <- new_data[order(new_data$num_id, new_data$int_time), ]
  # add between-level variables again
  new_data <- merge(
    x = new_data, y = btw_data, by = c("num_id"), all.x = TRUE, sort = FALSE
  )

  return(new_data)
}
