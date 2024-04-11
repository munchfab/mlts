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
#' @return A `data.frame` with missings imputed for use in `prepare_data()`.
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
  # create numeric id variable
  data$num_id <- as.numeric(as.factor(data[, id]))

  # remove between-level variables to add them later again
  btw_vars <- c(id, btw_vars)
  btw_data <- unique(data[, c("num_id", btw_vars)])
  data[, btw_vars] <- NULL

  for (i in 1:length(unique(data$num_id))) {

    # id
    num_id <- i

    # store continuous time variable for id == i
    time_cont <- data[data$num_id == i, time]

    # seq in steps of time grid from min to max of time
    time_seq <- seq(min(time_cont), max(time_cont), by = tinterval)

    # indicate shift interval
    shifts <- findInterval(x = time_cont, vec = time_seq)
    max_shifts <- max(table(shifts))

    # print warning if max_shifts increases 4
    if (max_shifts > 4) {
      warn_id <- btw_data[btw_data$num_id == i, id]
      warning("Observations shifted ", max_shifts, " intervals for id ", warn_id)
    }

    # create empty vector of imputed time
    time_imp <- NA

    # iterate through values of time_cont
    for (j in 1:length(time_cont)) {
      # fill time_imp with values of time_cont
      time_imp[
        # determine in which interval between two numbers of time_seq
        # lies value [j] of time_cont
        findInterval(x = time_cont[j], vec = time_seq)
      ] <- time_cont[j] # and fill that position with with value [j]
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
  new_data <- merge(x = imp_data, y = data, by = c("num_id", time), all = T)
  # sort by num_id and int_time
  new_data <- new_data[order(new_data$num_id, new_data$int_time), ]
  # add between-level variables again
  new_data <- merge(
    x = new_data, y = btw_data, by = c("num_id"), all.x = T, sort = F
  )

  return(new_data)
}


#' Create Missings for Approximation of Continuous Time Dynamic Models
#' (Deprecated)
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
#' @param clean Logical. Remove helper columns. Defaults to `TRUE`.
#'
#' @return A `data.frame` with missings imputed for use in `prepare_data()`.
#'
#' @export
create_missings2 <- function(data, tinterval, id, time,
                             btw_vars = NULL, clean = TRUE) {

  # create empty list for imputed data frames
  imp_list <- list()
  # create numeric id variable
  data$num_id <- as.numeric(as.factor(data[, id]))

  # remove between-level variables to add them later again
  btw_vars <- c(id, btw_vars)
  btw_data <- unique(data[, c("num_id", btw_vars)])
  data[, btw_vars] <- NULL

  # loop over observations in data set
  for (i in 1:length(unique(data$num_id))) {
    # store continuous time variable for id == i and have it start with 0
    diff_cont_time <- data[data$num_id == i, time] -
      min(data[data$num_id == i, time])
    # rescale to integer time variable as in Asparouhov et al. (2018)
    int_time <- ifelse(
      # if it rounds to 0 (for first occasion), make it 1 for first measurement
      ceiling(diff_cont_time / tinterval) == 0, 1,
      # otherwise just round up
      ceiling(diff_cont_time / tinterval)
    )
    temp <- data.frame(diff_cont_time, int_time)
    # calculate difference between integer time values, starting with 0
    temp$diff_int_time <- c(0, diff(int_time))
    # store at which position the difference of time_int is 0
    # (i.e., integer time is in same interval)
    overlap <- which(diff(int_time) == 0) + 1
    # store at which position the difference of time_int is greater 0
    # (i.e., integer time is not in same interval)
    nonoverlap <- which(diff(int_time) > 1) + 1
    # create indicator variables for (non-)overlapping time intervals
    temp[overlap, "overlap"] <- 1
    temp[nonoverlap, "nonoverlap"] <- 1
    temp$int_time_shifted <- temp$int_time
    # count number of times where observations are in same interval
    n_overlap <- length(overlap)
    # if there is overlap, shift the obversavations as in Hamaker (2018)
    if (n_overlap > 0) {
      # start counter
      k <- 1
      # repeat shifting until k > n_overlap
      repeat {
        min <- overlap[1]
        max <- nonoverlap[which(nonoverlap > overlap[1])][1]
        temp[min:(max - 1), "int_time_shifted"] <-
          temp[min:(max - 1), "int_time_shifted"] + 1
        overlap <- which(diff(temp$int_time_shifted) == 0) + 1
        nonoverlap <- which(diff(temp$int_time_shifted) > 1) + 1
        nonoverlap[length(nonoverlap) + 1] <- dim(temp)[1] + 1
        k <- k + 1
        if (k > n_overlap) {break}
      }
    }
    # check whether or not shifted more than 4 intervals and print warning
    temp$shifts <- temp$int_time_shifted - temp$int_time
    if (max(temp$shifts) > 4) {
      warning("Observations shifted", temp$shifts, "intervals for id", id)
    }
    # rescale to start counts with 1
    temp$int_time_shifted <- temp$int_time_shifted -
      min(temp$int_time_shifted) + 1
    temp[, time] <- data[data$num_id == i, time]
    # insert missings
    new <- data.frame(num_id = i,
                      int_time = seq(1:temp[dim(temp)[1], "int_time"]))

    # store in list
    imp_list[[i]] <- merge(new, temp, by = "int_time", all = TRUE)
  }
  # reduce to data frame
  imp_data <- do.call("rbind", imp_list)
  # join with existing data set
  # with base R, somewhat slower than with dplyr
  new_data <- merge(x = imp_data, y = data, by = c("num_id", time), all = T)
  # sort by num_id and int_time
  new_data <- new_data[order(new_data$num_id, new_data$int_time), ]
  # add between-level variables again
  new_data <- merge(
    x = new_data, y = btw_data, by = c("num_id"), all.x = T, sort = F
  )


  if (clean == TRUE) {
    # remove helper columns
    helpers <- c(
      "diff_cont_time", "diff_int_time", "overlap", "nonoverlap",
      "int_time_shifted", "shifts"
    )
    new_data <- new_data[, !colnames(new_data) %in% helpers]
  }

  return(new_data)
}

