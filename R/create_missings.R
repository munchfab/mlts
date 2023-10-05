create_data <- function(data, delta, id, time) {
  
  # create empty list for imputed data frames
  imp_list <- list()
  
  # create numeric id variable
  data$num_id <- as.numeric(as.factor(data[, id]))
  
  # loop over obersavtions in data set
  for (i in 1:length(unique(data$num_id))) {
    # store continuous time variable for id == i and have it start with 0
    diff_cont_time <- data[data$num_id == i, time] -
      min(data[data$num_id == i, time])
    # rescale to integer time variable as in Hamaker et al. (2018)
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
    # imp_list[[i]] <- dplyr::full_join(new, temp, by = "int_time")
    # imp_list[[i]] <- temp
  }
  
  # reduce to data frame
  imp_data <- purrr::reduce(imp_list, rbind)
  
  # join with existing data set
  stan_data <- dplyr::full_join(
    imp_data, data,
    by = c("num_id", time)
  )
  # fill original IDs up
  stan_data <- dplyr::group_by(.data = stan_data, num_id)
  stan_data <- tidyr::fill(stan_data, all_of(id), .direction = "downup")
  stan_data <- dplyr::ungroup(stan_data)
  
  return(stan_data)
}
