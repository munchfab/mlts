#' Prepare data for Stan
#'
#' @param data An object of class `data.frame` (or one that can be coerced
#' to that class) in long format.
#' @param id Character. The variable in `data` that identifies the person or observational
#' unit.
#' @param ts Character. The variable(s) in `data` that
#' contain the time-series construct. If multiple variables are provided in a
#' character vector, a vector autoregressive model is fit.
#' @param time Character. The variable in `data` that contains the (continuous) time.
#' @param tinterval The step interval for approximation for a continuous time
#' dynamic model. The smaller the step interval, the better the approximation.
#' @param beep Character. The variable in `data` that contains the running
#' beep from 1 to TP for each person.
#' @param days Optional. If a running beep identifier is provided via the `beep`
#' argument and observations are nested within days (or similar grouping unit),
#' the variable in `data` that contains the day identifier can be added to correct
#' for overnight lags (see Details).
#' @param n_overnight_NAs Optional. The number of `NA` rows to add after the last
#' observation of each day (if `day` is provided).
#' @param na.rm logical. As default option missing values remain in the data and
#' will be imputed during model estimation. Set to `TRUE` to remove all rows with
#' missing values in variables given in `ts`.
#' @param covariates Named character vector. An optional named vector of
#' characters to refer to predictors of random effects as specified in the `model`.
#' Note that specifying `covariates` is only necessary if the respective
#' variable name(s) in `data` differ from the variables names specified in `model`.
#' @param outcomes Named character vector. Similar to `covariates`, an optional named vector of
#' characters to refer to outcome predicted by random effects as specified in the `model`.
#' Note that specifying `outcomes` is only necessary if the respective
#' variable name(s) in `data` differ from the outcome variable name(s) specified in `model`.
#' @param outcome_pred_btw character.
#'
#' @return A `data.frame` that can be passed to \code{\link[mlts]{mlts_fit}}.
#' @noRd
#'
#' @examples
#' \donttest{
#' # prepare data for vector-autoregressive model
#' data <- prepare_data(
#'   ts = c("Y1", "Y2"),
#'   data = ts_data,
#'   id = "ID",
#'   time = "time",
#'   tinterval = 1,
#' )
#'
#' # further examples with overnight lags
#'
#' }
prepare_data <- function(data, id, ts, time = NULL, tinterval = NULL,
                         beep = NULL, days = NULL,
                         n_overnight_NAs, na.rm = FALSE, covariates = NULL,
                         outcomes = NULL, outcome_pred_btw = NULL){

  # coerce to data frame if necessary (e.g., if tibble is provided)
  data <- as.data.frame(data)

  # create a subset of the data with fixed variable names
  data$num_id = data[,id]

  # ensure that the id-variable is running from 1 to N
  data$num_id = as.numeric(factor(data$num_id))

  # depending on inputs
  if(!is.null(time)){
    data$time = data[,time]
    data = data[order(data$num_id, data$time),]
  }
  if(!is.null(beep)){
    data$beep = data[,beep]
    data = data[order(data$num_id, data$beep),]
  }
  if(!is.null(days)){
    data$day = data[,days]
    data = data[order(data$num_id, data$day, data$beep),]
  }

  # store between person variables
  btw.vars = c(names(outcomes), names(covariates), names(outcome_pred_btw))

  # general step: add variable to detect observations with NAs
  # remove rows containing missing values
  if(length(ts) == 1){
    data$miss.NA = is.na(data[,ts])
  } else {
    data$miss.NA = rowSums(is.na(data[,ts])) > 0
  }

  # General first step: ======================================================
  # remove all (consecutive) NAs on first or last observations
  N_obs_id = data.frame(table(data$num_id))$Freq   # number of obs per subject
  data$obs_number <- NA
  data$obs_number_rev <- NA
  data$miss.NA.cumsum <- NA
  data$miss.NA.cumsum_rev <- NA
  for(i in 1:max(data$num_id)){
    data$obs_number[data$num_id == i] = 1:N_obs_id[i]
    data$miss.NA.cumsum[data$num_id == i] = cumsum(data$miss.NA[data$num_id == i])
    data$obs_number_rev[data$num_id == i] = N_obs_id[i]:1
    data$miss.NA.cumsum_rev[data$num_id == i] = rev(cumsum(rev(data$miss.NA[data$num_id == i])))
  }

    # mark rows to exclude
    data$excl = ifelse((data$obs_number == data$miss.NA.cumsum) |
                       (data$obs_number_rev == data$miss.NA.cumsum_rev), 1, 0)
    data = data[data$excl == 0,]

    # remove helper columns
    data[,c("obs_number", "obs_number_rev",
            "miss.NA.cumsum", "miss.NA.cumsum_rev", "excl")] <- NULL


  # ===========================================================================

  # depending on inputs
  if(na.rm == TRUE){
    # remove rows containing missing values
    data <- data[data$miss.NA == FALSE,]

    # print warning for NA removal and tinterval
    if (is.numeric(tinterval)) {
      warning("Removing NAs with \"na.rm = TRUE\" does not allow to approximate",
              " a continuous time model.")
    }

  } else if (is.numeric(tinterval)) {
    if (is.character(time)) {
      # create time grid according to continuous time variable
      data = create_missings(data = data, tinterval = tinterval, id = id,
                             time = time, btw_vars = btw.vars)
      data = cbind("order" = 1:nrow(data), data)

      # ADD LATER ==============================================================
      # check number of missings and print a warning in case number is high

      # ========================================================================

      # merge data to time-grid to avoid missing values in between-person variables
      data = data[order(data$order),]

      # to speed up estimation: Exclude all (consecutive) NAs on first observations
    } else if (is.character(days)) {
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
      data = create_missings(data = data, tinterval = tinterval, id = id,
                             time = "beep_new", btw_vars = btw.vars)
      data = cbind("order" = 1:nrow(data),data)

      # ADD LATER ==============================================================
      # check number of missings and print a warning in case number is high

      # ========================================================================
      data = data[order(data$order),]
    } else {
      ############################ CONTINUE HERE FOR STANDARD BEEP OPTION -----

      # stop and print warning for missing time or beep variable
      stop("Specifying a tinterval requires additional specification of",
           " either a \"time\" or a \"beep\" variable in data.")
    }
  }


    # remove helper column(s)
    data[,c("miss.NA")] <- NULL


  return(data)
}
