#' Create data for DSEM models
#'
#' @param data An object of class `data.frame` (or one that can be coerced
#' to that class) containing data of all variables used in the model.
#' @param delta The step interval for approximation for a continuous time DSEM.
#' The smaller the step interval, the better the approximation.
#' @param id The variable in `data` that identifies the person or observational
#' unit (as string).
#' @param time The variable in `data` that contains the (continuous) time (as
#' string).
#'
#' @return A `data.frame` with missings imputed for use in `dsem()`.
#'
#' @examples Tbd
#' @export
create_dsem_data <- function(data, delta, id, time) {

  # create empty list for imputed data frames
  imp_list <- list()

  # create numeric id variable
  data$num_id <- as.numeric(as.factor(data[, id]))

  # loop over obersavtions in data set
  for (i in 1:length(unique(data$num_id))) {
    # store continuous time variable for id == i and have it start with 0
    diff_cont_time <- data[data$num_id == i, time] -
      min(data[data$num_id == i, time])
    # rescale to integer time variable as in Asparouhov et al. (2018)
    int_time <- ifelse(
      # if it rounds to 0 (for first occasion), make it 1 for first measurement
      ceiling(diff_cont_time / delta) == 0, 1,
      # otherwise just round up
      ceiling(diff_cont_time / delta)
    )
    # create temporary data frame continuous and integer time variables
    temp <- data.frame(diff_cont_time, int_time)
    # calculate difference between integer time values, starting with 0
    temp$diff_int_time <- c(0, diff(int_time))
    # store at which position the difference of time_int is 0
    # (i.e., integer time is in same interval)
    overlap <- which(diff(int_time) == 0) + 1
    # store at which position the difference of time_int is greater 0
    # (i.e., integer time is not in same interval)
    nonoverlap <- which(diff(int_time) > 1) + 1
    # create indicator variables for overlapping time intervals
    temp[overlap, "overlap"] <- 1
    # create indicator variables for non-overlapping time intervals
    temp[nonoverlap, "nonoverlap"] <- 1
    # create new shifted time variable fo counting number of shifts
    temp$int_time_shifted <- temp$int_time
    # count number of times where observations are in same interval
    n_overlap <- length(overlap)

    # if there is overlap, shift the obversavations as in Hamaker (2018)
    if (n_overlap > 0) {
      # start counter
      k <- 1
      # repeat shifting until k > n_overlap
      repeat {
        # no idea whats happening here
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
    # check whether or not shifted more than 4 intervals
    temp$shifts <- temp$int_time_shifted - temp$int_time
    # if so, print warning
    if (max(temp$shifts) > 4) {
      print(paste("Warning: obs shifted", temp$shifts,
                  "intervals for id", id))
    }
    # rescale to start counts with 1
    temp$int_time_shifted <- temp$int_time_shifted -
      min(temp$int_time_shifted) + 1
    temp[, time] <- data[data$num_id == i, time]
    # insert missing
    new <- data.frame(num_id = i,
                      int_time = seq(1:temp[dim(temp)[1], "int_time"]))

    # store in list
    imp_list[[i]] <- merge(new, temp, by = "int_time", all = TRUE)
  }

  # reduce to data frame
  imp_data <- purrr::reduce(imp_list, rbind)

  # join with existing data set
  dsem_data <- dplyr::full_join(
    imp_data, data,
    by = c("num_id", time)
  )
  # fill original IDs up
  sem_data <- dplyr::group_by(.data = dsem_data, num_id)
  dsem_data <- tidyr::fill(dsem_data, id, .direction = "downup")
  dsem_data <- dplyr::ungroup(dsem_data)
  # want to do with pipe but doesnt work
  # dsem_data <- dplyr::group_by(.data = dsem_data, num_id) %>%
  #   tidyr::fill(id, .direction = "downup") %>%
  #   dplyr::ungroup(dsem_data)

  return(dsem_data)
}
