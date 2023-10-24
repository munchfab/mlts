#' Create data for Stan
#'
#' @param data An object of class `data.frame` (or one that can be coerced
#' to that class) in long format.
#' @param id character. The variable in `data` that identifies the person or
#' observational unit.
#' @param beep character. The variable in `data` that contains the running
#' beep from 1 to TP for each person.
#' @param ts character. The variable in `data` that contains the observations
#' of the time-varying construct.
#' @param pred_random character. Name(s) of between-person variables to use as
#' predictors of individual parameters
#' @param outcome character. Names(s) of between-person variables
#' to regress individual parameters on.
#' @param out_predictors character. "mu", "ar", and/or "logv" to use individual
#' parameter estimates.
#' Can be combined with any other between-person variable.
#' @param standardize_out logical. Should output be standardized? Defaults
#' to `TRUE`.
#' @param overnight_lags character. Add column name of the day number to avoid
#' using the last observation of a day as lagged predictor of the
#' subsequent observation on the next day.
#' @param miss_handling character. Should missings be removed (`"remove"`) or
#' imputed (`"impute"`)?
#' @param random.innovations logical. `FALSE` = constant innovation variance,
#' `TRUE` = person-specific innovation variances.
#' @param add_mplus_data logical. Should data be stored in Mplus format?
#' Defaults to `TRUE`.
#'
#' @return A `list` that can be passed to `stan()`.
#' @export
#'
#' @examples TBA
create_stan_data <- function(data, id, beep, ts,
                             pred_random = NULL,
                             outcome = NULL,
                             out_predictors = NULL,
                             standardize_out = TRUE,
                             overnight_lags = NULL,
                             miss_handling = c("remove", "impute"),
                             random.innovations = TRUE,
                             add_mplus_data = T) {


  # integrate later...
  out_pred_b = NULL
  # create copy of data with unique names
  df <- data.frame("id" = data[, id],
                   "dayno" = 1,
                   "beep" = data[, beep])

  # add observations
  df = cbind(df, data[, ts])
  # add day identifier if provided
  if (!is.null(overnight_lags)) {
    df$dayno = data[, overnight_lags]
  }

  # add between-person variables if provided
  if (!is.null(pred_random) | !is.null(outcome)) {
    out_pred_b = out_predictors[!(out_predictors %in% c("mu", "ar", "logv"))]
    if (length(out_pred_b) == 0) {
      df = cbind(df, data[, c(pred_random, outcome)])
      out_pred_b = NULL
    } else {
      df = cbind(df, data[, c(pred_random, outcome, out_pred_b)])
    }
  }

  # make sure that old variable names are removed
  # in case old column names were passed to df
  if (length(ts) > 1) {ts_col = paste0("ts_", 1:length(ts))}
  else {ts_col = "ts"}
  colnames(df) = c(
    "id", "dayno", "beep",
    ts_col, pred_random, outcome, out_pred_b
  )

  # handling of missing values
  if (miss_handling == "remove") {
    # remove all missings
    df <- dplyr::filter(df, dplyr::if_all(dplyr::everything(), ~ !is.na(.)))
  } else if (miss.handling == "impute") {
    # set all missings to -99
    df <- dplyr::mutate(
      df,
      dplyr::across(.cols = grepl(pattern = "ts", colnames(.)),
                    ~ ifelse(is.na(.), -99, .))
    )
  }
}
