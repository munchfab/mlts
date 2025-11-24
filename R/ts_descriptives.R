#' Describe Time Series Variables by Group
#'
#' Computes descriptive statistics for one or more time-series variables,
#' optionally within groups. For each variable (and group, if provided), the
#' function returns sample size, number of complete observations, percentage of
#' missing values, number of complete lag-1 pairs, central tendency, dispersion,
#' skewness, kurtosis, and floor/ceiling effects.
#'
#' @param data A `data.frame` containing the time-series variables.
#' @param ts A character vector of variable names (columns in \code{data})
#'   representing the time-series to be summarized.
#' @param group Optional. A character string giving the name of a grouping
#'   variable in \code{data}. If \code{NULL} (default), all observations are
#'   treated as belonging to one group.
#' @param digits Integer indicating the number of decimal places to which numeric
#'   statistics should be rounded. Default is 3.
#' @param scale_min,scale_max Optional numeric values specifying the minimum and
#'   maximum of the measurement scale for computing floor and ceiling
#'   percentages. If \code{NULL} (default), these are taken as the minimum and
#*   maximum observed values in the entire dataset for each variable.
#'
#' @details
#' For each time-series variable, the function removes \code{NA} values before
#' computing summary statistics. Skewness and kurtosis are computed manually and
#' set to \code{NA} when the standard deviation is zero, missing, or when fewer
#' than four complete observations are available.
#'
#' Floor and ceiling effects are defined as the percentage of non-missing
#' observations equal to \code{scale_min} and \code{scale_max}, respectively.
#'
#' Lag-1 completeness is computed as the number of adjacent observation pairs
#' where both the current and next observation are non-missing. **If
#' \code{group = NULL}, this computation does not account for clustering
#' or overnight gaps**—that is, the function assumes the
#' data form a single continuous time series, and lag-1 pairs may span across
#' natural breaks in the data.
#'
#' @return
#' A `data.frame` with one row per variable per group, containing:
#' \itemize{
#'   \item \code{group} — group identifier (omitted if \code{group = NULL})
#'   \item \code{ts_var} — name of time series variable
#'   \item \code{N} — total number of observations (including NAs)
#'   \item \code{N_comp} — number of non-missing observations
#'   \item \code{miss_pc} — percent missing
#'   \item \code{N_lag1_comp} — number of complete lag-1 pairs
#'   \item \code{M} — mean
#'   \item \code{Mdn} — median
#'   \item \code{SD} — standard deviation
#'   \item \code{Skew} — skewness
#'   \item \code{Kurtosis} — kurtosis (excess)
#'   \item \code{floor_pc} — percent at \code{scale_min}
#'   \item \code{ceiling_pc} — percent at \code{scale_max}
#' }
#'
#' @examples
#' \dontrun{
#' # Example with grouping
#' describe_ts(data = df, ts = c("x", "y"), group = "id")
#'
#' # Example without grouping
#' describe_ts(data = df, ts = "x")
#' }
#'
#' @export
describe_ts <- function(data, ts, group = NULL, digits = 3, scale_min = NULL, scale_max = NULL){

  if(!is.null(group)){
    data$temp_id <- data[[group]]
  } else {
    data$temp_id <- 1              # pseudo-id for grouping
  }
  groups <- unique(data$temp_id)


  results <- list()
  row_index <- 1

  for (gg in groups){
    data_p <- data[ data$temp_id == gg, ]

    for (var in ts){

      if(is.null(scale_min)){
        scale_min = min(data_p[[var]], na.rm = T)
      }
      if(is.null(scale_max)){
        scale_max = max(data_p[[var]], na.rm = T)
      }

      x_vec <- data_p[[var]]
      lag_pos = 1:(length(x_vec)-1)
      n_lag_complete = sum(!is.na(x_vec[lag_pos+1]) & !is.na(x_vec[lag_pos]))

      na_missing_num <- sum(is.na(x_vec))
      n_with_na    <- length(x_vec)
      na_missing_p <- (na_missing_num / n_with_na) * 100

      x_vec_clean <- x_vec[!is.na(x_vec)]
      n           <- length(x_vec_clean)
      x_mean      <- mean(x_vec_clean)
      x_median    <- stats::median(x_vec_clean)
      x_sd        <- stats::sd(x_vec_clean)
      pc_floor    <- sum(x_vec_clean == scale_min) / n * 100
      pc_ceiling  <- sum(x_vec_clean == scale_max) / n * 100



      if (x_sd == 0 || is.na(x_sd) || n < 4){
        x_kurt <- NA_real_
        x_skew <- NA_real_
      } else{
        x_kurt <- sum(((x_vec_clean - x_mean) / x_sd)^4) / n - 3
        x_skew <- sum(((x_vec_clean - x_mean) / x_sd)^3) / n
      }

      one_row <- list(
        group      = gg,
        ts_var     = var,
        N          = n_with_na,
        N_comp     = n,
        miss_pc    = na_missing_p,
        N_lag1_comp= n_lag_complete,
        M          = x_mean,
        Mdn        = x_median,
        SD         = x_sd,
        Skew       = x_skew,
        Kurtosis   = x_kurt,
        floor_pc   = pc_floor,
        ceiling_pc = pc_ceiling
      )
      results[[row_index]] <- one_row
      row_index <- row_index + 1
    }
  }

  final_df <- do.call(rbind, lapply(results, as.data.frame))
  final_df[5:ncol(final_df)] <- apply(final_df[5:ncol(final_df)], 2, function(x) round(x, digits))


  if(is.null(group)){
    final_df$group <- NULL
  }


  return(final_df)

}
