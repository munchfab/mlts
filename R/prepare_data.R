#' Prepare data for Stan
#'
#' @param data An object of class `data.frame` (or one that can be coerced
#' to that class) in long format.
#' @param id character. The variable in `data` that identifies the person or
#' observational unit.
#' @param ts_ind character. The variable(s) in `data` that contain(s) the time series data.
#' @param time character. The variable in `data` that contains the (continuous) time (as
#' string).
#' @param tinterval The step interval for approximation for a continuous time DSEM.
#' The smaller the step interval, the better the approximation.
#' @param beep character. The variable in `data` that contains the running
#' beep from 1 to TP for each person.
#' @param days Optional. If a running beep identifier is provided via the `beep` argument and
#' observations are nested within days (or similar grouping unit), the variable
#' in `data` that contains the day identifier can be added to correct for overnight lags (see Details).
#' @param n_overnight.NAs Optional. The number of `NA` rows to add after the last observation of each day (if `day` is provided).
#' @param na.rm logical. As default option missing values remain in the data and
#' will be imputed during model estimation. Set to `TRUE` to remove all rows with
#' missing values in variables given in `ts_ind`.
#'
#' @return A `list` that can be passed to `stan()`.
#' @export
#'
#' @examples TBA
prepare_data <- function(data, id, ts_ind, time, tinterval, beep, days,
                         n_overnight_NAs, na.rm = FALSE, outcomes = NULL,
                         covariates = NULL, outcome.pred.btw = NULL){


  # create a subset of the data with fixed variable names
  data$num_id = data[,id]

  # depending on inputs
  if(is.character(time)){
    data$time = data[,time]
    data = data[order(data$num_id, data$time),]
  }
  if(is.character(beep)){
    data$beep = data[,beep]
    data = data[order(data$num_id, data$beep),]
  }
  if(is.character(days)){
    data$day = data[,days]
    data = data[order(data$num_id, data$day, data$beep),]
  }

  # ensure that the id-variable is running from 1 to N
  data$num_id = as.numeric(factor(data$num_id))

  # store between person variables
  btw.vars = c(names(outcomes), names(covariates), names(outcome.pred.btw))

  # depending on inputs
  if(na.rm == TRUE){
    # remove rows containing missing values
    data = data[rowSums(is.na(data[,ts_ind])) == 0, ]

    } else if(is.character(time) & is.numeric(tinterval)){
      # create time grid according to continuous time variable
      data = create_missings(data = data, delta = tinterval, id = id,
                      time = "time", btw.vars = btw.vars)
      data = cbind("order" = 1:nrow(data), data)

      # ADD LATER ==============================================================
      # check number of missings and print a warning in case number is high

      # ========================================================================

      # merge data to time-grid to avoid missing values in between-person variables
      data = data[order(data$order),]

      } else if(is.character(days) & is.numeric(tinterval)){
      # create a pseudo time-grid using the beep number for data nested within days
        if(!is.character(beep)){
          stop("No beep variable provided")
        }
        if(tinterval != 1){
          warning("tinterval should be 1")
        }
        if(!is.numeric(n_overnight_NAs)){
          stop("Specify the number of missings to add.")
        }
        # update the beep numbers by a multiple of the day number and n_overnight_NAs
        add_beep = (data$day-1) * n_overnight_NAs
        new_beep = add_beep + data$beep
        data$beep_new = new_beep

        # create time grid according to continuous time variable
        data = create_missings(data = data, delta = tinterval, id = id,
                                    time = "beep_new", btw.vars = btw.vars)
        data = cbind("order" = 1:nrow(data),data)

        # ADD LATER ==============================================================
        # check number of missings and print a warning in case number is high

        # ========================================================================
        data = data[order(data$order),]

    } else {
      ############################ CONTINUE HERE FOR STANDARD BEEP OPTION --------

    }

  return(data)
}
