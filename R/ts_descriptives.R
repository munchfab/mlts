


AR1_4ind = mlts_model(
  q = 1,          # number of time-varying (latent) constructs
  p = 4,          # number of manifest indicators
  btw_factor = T  # common between-level factor (the default)
)

simAR1_4ind = mlts_sim(
  model = AR1_4ind,         # model object
  default = T,              # set to TRUE to use random true values for parameters
  N = 100,                  # number of subjects
  TP = 70,                  # number of time points per subject
  burn.in = 50,             # burn in for within-level process
  seed = 4,                 # set seed for reproducible results
)

# testing
df_test <- data.frame(
  a = c(1, 2, 3, 4, 5),
  b = c(10, 10, 10, 10, 10),
  c = c(NA, 1, 1, 1, 1)
)



# Input(data frame with data, vector with the ts - variables)
ts_descriptives <- function(data, ts, per_person = F, id = NULL){

  if (per_person == F){
    # number of observations
    n_obs = sapply(ts, function(x) sum(!is.na(data[[x]])))

    # mean
    mean_ts = sapply(ts, function(x) mean(data[[x]], na.rm = T))

    # median
    median_ts = sapply(ts, function(x) stats::median(data[[x]], na.rm = T))

    # sample standard deviation
    sd_ts = sapply(ts, function(x) {
      x_vec <- data[[x]]
      n <- sum(!is.na(x_vec))

      x_var <- sum((x_vec - mean(x_vec, na.rm = T))^2, na.rm = T) / (n - 1)
      sqrt(x_var)
    })

    # kurtosis
    kurt_ts <- sapply(ts, function(x) {
      x_vec <- data[[x]]
      x_vec <- x_vec[!is.na(x_vec)]
      n <- length(x_vec)

      x_mean <- mean(x_vec)
      x_var <- sum((x_vec - x_mean)^2) / (n - 1)
      x_sd <- sqrt(x_var)

      if (n < 4 || x_sd == 0) return(NA_real_)

      # Moment Kurtosis
      x_kurt <- (sum(((x_vec - x_mean) / x_sd)^4)) / n

      # Excess Kurtosis
      x_kurt - 3
    })

    # Moment Skewness

    skew_ts <- sapply(ts, function(x) {
      x_vec <- data[[x]]
      x_vec <- x_vec[!is.na(x_vec)]
      n <- length(x_vec)

      x_mean <- mean(x_vec)
      x_var <- sum((x_vec - x_mean)^2) / (n - 1)
      x_sd <- sqrt(x_var)

      if(n < 3 || x_sd == 0) return(NA_real_)

      x_skew <- sum(((x_vec - x_mean) / x_sd)^3) / n
    })

    # fehlende Werte
    n_missing_num <- sapply(ts, function(x) sum(is.na(data[[x]])))
    n_with_na <- sapply(ts, function(x) length(data[[x]]))
    n_missing_p <- sapply(ts, function(x) {
      (sum(is.na(data[[x]])) / length(data[[x]])) * 100
    })

    df <- data.frame(variable = ts, n_obs = n_obs, mean = mean_ts,
                     median = median_ts, sd = sd_ts,
                     skew = skew_ts, kurt = kurt_ts,
                     missing_n = n_missing_num, missing_p = n_missing_p)
    rownames(df) <- NULL
    return(df)
  }


  else {
    persons <- unique(data[[id]]) # Personen
    results <- list()
    row_index <- 1

    for (person in persons){ # looping through every person
      data_p <- data[data[[id]] == person, ] # list of each persons data

      for (var in ts){ # looping through every variable for each person
        x_vec <- data_p[[var]] # creating data for each variable

        # number of nas
        na_missing_num <- sum(is.na(x_vec))
        n_with_na <- length(x_vec)
        if (n_with_na == 0){ # check for 0 values in variable x
          na_missing_p <- NA
        }
        else{
          na_missing_p <- (na_missing_num / n_with_na) *100
        }

        # nas entfernen
        x_vec_clean <- x_vec[!is.na(x_vec)]

        # number of observations
        n <- length(x_vec_clean)

        #mean
        x_mean <- mean(x_vec_clean)

        # median
        x_median <- stats::median(x_vec_clean)

        # sample stadard devation
        if (n < 2){
          x_sd <- NA
        }
        else {
          x_var <- sum((x_vec_clean - x_mean)^2) / (n - 1)
          x_sd <- sqrt(x_var)
        }
        if (x_sd == 0 || is.na(x_sd) || n < 4){ # consevative for skew (normally: n < 3)
          x_kurt <- NA
          x_skew <- NA
        }


        else{
          # Excess kurtosis
          x_kurt <- sum(((x_vec_clean - x_mean) / x_sd)^4) / n - 3
          #skew
          x_skew <- sum(((x_vec_clean - x_mean) / x_sd)^3) / n
        }
        one_row <- list(ID = person, variable = var,
                n_obs = n, mean = x_mean, median = x_median,
                sd = x_sd, skewness = x_skew, kurtosis = x_kurt,
                missing_n = na_missing_num, missing_p = na_missing_p)
        results[[row_index]] <- one_row
        row_index <- row_index + 1
      }
    }
    final_df <- do.call(rbind, lapply(results, as.data.frame))
    return(final_df)
  }

}

ts_descriptives(simAR1_4ind$data, c("Y1.1", "Y1.2", "Y1.3", "Y1.4"))

ts_descriptives(simAR1_4ind$data, c("Y1.1", "Y1.2", "Y1.3", "Y1.4"),
                per_person = T, id = "ID")
ts_descriptives(df_test, c("a", "b", "c"))
