#' Descriptive statistics for time-series variables
#'
#' Computes univariate descriptive statistics (mean, median, standard deviation,
#' skewness, excess kurtosis, and amount of missing data) for one or more
#' time-series variables, optionally separately for each individual.
#'
#' @param data A \code{data.frame} containing the time-series data.
#' @param ts A character vector with the names of the time-series variables.
#' @param per_person Logical. If \code{FALSE} (default), statistics are
#'   computed across all available observations in \code{data}. If
#'   \code{TRUE}, statistics are computed separately for each individual
#'   identified by \code{id}.
#' @param id Optional character string giving the name of the ID variable in
#'   \code{data}. Required if \code{per_person = TRUE}.
#'
#' @details
#' For each variable in \code{ts}, the function uses all non-missing
#' observations to compute:
#' \itemize{
#' \item \code{n_obs}: number of non-missing observations
#' \item \code{mean}: arithmetic mean
#' \item \code{median}: sample median
#' \item \code{sd}: sample standard deviation (with denominator \code{n - 1})
#' \item \code{skew}/\code{skewness}: moment-based skewness
#' \item \code{kurt}/\code{kurtosis}: excess kurtosis (moment-based, i.e. kurtosis - 3)
#' \item \code{missing_n}: number of missing values
#' \item \code{missing_p}: percentage of missing values
#' }
#'
#' For small sample sizes or zero variance, skewness and kurtosis are returned
#' as \code{NA}. Specifically, for the overall descriptives,
#' kurtosis is set to \code{NA} if \code{n < 4} or the standard deviation is
#' zero; skewness is set to \code{NA} if \code{n < 3} or the standard deviation
#' is zero. The same rules are applied in the per-person case.
#'
#' @return
#' If \code{per_person = FALSE}, a \code{data.frame} with one row per variable
#' in \code{ts} and columns
#' \code{variable}, \code{n_obs}, \code{mean}, \code{median}, \code{sd},
#' \code{skew}, \code{kurt}, \code{missing_n}, \code{missing_p}.
#'
#' If \code{per_person = TRUE}, a \code{data.frame} with one row per
#' person-variable combination and columns
#' \code{ID}, \code{variable}, \code{n_obs}, \code{mean}, \code{median},
#' \code{sd}, \code{skewness}, \code{kurtosis}, \code{missing_n},
#' \code{missing_p}.
#'
#' @examples
#' ## Example data
#' set.seed(1)
#' df <- data.frame(
#'   ID = rep(1:2, each = 5),
#'   y1 = rnorm(10),
#'   y2 = rnorm(10)
#' )
#'
#' ## Overall descriptives across all persons
#' ts_descriptives(df, ts = c("y1", "y2"))
#'
#' ## Per-person descriptives
#' ts_descriptives(df, ts = c("y1", "y2"), per_person = TRUE, id = "ID")
#'
#' @export
ts_descriptives <- function(data, ts, per_person = FALSE, id = NULL){

  if (per_person == FALSE){
    # number of observations
    n_obs = sapply(ts, function(x) sum(!is.na(data[[x]])))

    # mean
    mean_ts = sapply(ts, function(x) mean(data[[x]], na.rm = TRUE))

    # median
    median_ts = sapply(ts, function(x) stats::median(data[[x]], na.rm = TRUE))

    # sample standard deviation
    sd_ts = sapply(ts, function(x) {
      x_vec <- data[[x]]
      n <- sum(!is.na(x_vec))

      x_var <- sum((x_vec - mean(x_vec, na.rm = TRUE))^2, na.rm = TRUE) / (n - 1)
      sqrt(x_var)
    })

    # kurtosis (excess)
    kurt_ts <- sapply(ts, function(x) {
      x_vec <- data[[x]]
      x_vec <- x_vec[!is.na(x_vec)]
      n <- length(x_vec)

      x_mean <- mean(x_vec)
      x_var <- sum((x_vec - x_mean)^2) / (n - 1)
      x_sd <- sqrt(x_var)

      if (n < 4 || x_sd == 0) return(NA_real_)

      # Moment kurtosis
      x_kurt <- sum(((x_vec - x_mean) / x_sd)^4) / n

      # Excess kurtosis
      x_kurt - 3
    })

    # skewness (moment-based)
    skew_ts <- sapply(ts, function(x) {
      x_vec <- data[[x]]
      x_vec <- x_vec[!is.na(x_vec)]
      n <- length(x_vec)

      x_mean <- mean(x_vec)
      x_var <- sum((x_vec - x_mean)^2) / (n - 1)
      x_sd <- sqrt(x_var)

      if (n < 3 || x_sd == 0) return(NA_real_)

      sum(((x_vec - x_mean) / x_sd)^3) / n
    })

    # missing values
    n_missing_num <- sapply(ts, function(x) sum(is.na(data[[x]])))
    n_missing_p <- sapply(ts, function(x) {
      (sum(is.na(data[[x]])) / length(data[[x]])) * 100
    })

    df <- data.frame(
      variable   = ts,
      n_obs      = n_obs,
      mean       = mean_ts,
      median     = median_ts,
      sd         = sd_ts,
      skew       = skew_ts,
      kurt       = kurt_ts,
      missing_n  = n_missing_num,
      missing_p  = n_missing_p
    )
    rownames(df) <- NULL
    return(df)
  } else {
    persons <- unique(data[[id]]) # Personen
    results <- list()
    row_index <- 1

    for (person in persons){ # looping through every person
      data_p <- data[data[[id]] == person, ] # data for this person

      for (var in ts){ # looping through every variable for each person
        x_vec <- data_p[[var]] # data for each variable

        # number of NAs
        na_missing_num <- sum(is.na(x_vec))
        n_with_na <- length(x_vec)
        if (n_with_na == 0){
          na_missing_p <- NA_real_
        } else{
          na_missing_p <- (na_missing_num / n_with_na) * 100
        }

        # remove NAs
        x_vec_clean <- x_vec[!is.na(x_vec)]

        # number of observations
        n <- length(x_vec_clean)

        # mean
        x_mean <- mean(x_vec_clean)

        # median
        x_median <- stats::median(x_vec_clean)

        # sample standard deviation
        if (n < 2){
          x_sd <- NA_real_
        } else {
          x_var <- sum((x_vec_clean - x_mean)^2) / (n - 1)
          x_sd <- sqrt(x_var)
        }

        if (x_sd == 0 || is.na(x_sd) || n < 4){
          x_kurt <- NA_real_
          x_skew <- NA_real_
        } else{
          # excess kurtosis
          x_kurt <- sum(((x_vec_clean - x_mean) / x_sd)^4) / n - 3
          # skewness
          x_skew <- sum(((x_vec_clean - x_mean) / x_sd)^3) / n
        }

        one_row <- list(
          ID         = person,
          variable   = var,
          n_obs      = n,
          mean       = x_mean,
          median     = x_median,
          sd         = x_sd,
          skewness   = x_skew,
          kurtosis   = x_kurt,
          missing_n  = na_missing_num,
          missing_p  = na_missing_p
        )
        results[[row_index]] <- one_row
        row_index <- row_index + 1
      }
    }
    final_df <- do.call(rbind, lapply(results, as.data.frame))
    return(final_df)
  }
}


# AR1_4ind = mlts_model(
#   q = 1,          # number of time-varying (latent) constructs
#   p = 4,          # number of manifest indicators
#   btw_factor = T  # common between-level factor (the default)
# )
#
# simAR1_4ind = mlts_sim(
#   model = AR1_4ind,         # model object
#   default = T,              # set to TRUE to use random true values for parameters
#   N = 100,                  # number of subjects
#  TP = 70,                  # number of time points per subject
#   burn.in = 50,             # burn in for within-level process
#   seed = 4,                 # set seed for reproducible results
# )
#
#
# df_test <- data.frame(
#   a = c(1, 2, 3, 4, 5),
#   b = c(10, 10, 10, 10, 10),
#   c = c(NA, 1, 1, 1, 1)
# )
#
# ts_descriptives(simAR1_4ind$data, c("Y1.1", "Y1.2", "Y1.3", "Y1.4"))
#
# ts_descriptives(simAR1_4ind$data, c("Y1.1", "Y1.2", "Y1.3", "Y1.4"),
#                 per_person = T, id = "ID")
# ts_descriptives(df_test, c("a", "b", "c"))
