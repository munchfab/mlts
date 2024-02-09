#' Create data for Stan
#'
#' @param data An object of class `data.frame` (or one that can be coerced
#' to that class) in long format.
#' @param id character. The variable in `data` that identifies the person or
#' observational unit.
#' @param beep character. The variable in `data` that contains the running
#' beep from 1 to TP for each person.
#' @param y character. The variable in `data` that contains the observations
#' of the time-varying construct.
#' @param pred_random character. Name(s) of between-person variables to use as
#' predictors of individual parameters
#' @param outcome character. Name(s) of between-person variables
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
#' @param random_innovations logical. `FALSE` = constant innovation variance,
#' `TRUE` = person-specific innovation variances.
#' @param add_mplus_data logical. Should data be stored in Mplus format?
#' Defaults to `TRUE`.
#'
#' @return A `list` that can be passed to `stan()`.
#' @export
#'
#' @examples TBA
create_stan_data <- function(data, id, beep, y,
                             pred_random = NULL,
                             outcome = NULL,
                             out_predictors = NULL,
                             standardize_out = TRUE,
                             overnight_lags = NULL,
                             miss_handling = c("remove", "impute"),
                             random_innovations = TRUE,
                             add_mplus_data = TRUE) {


  # integrate later...
  out_pred_b = NULL
  # create copy of data with unique names
  df <- data.frame("id" = data[, id],
                   "dayno" = 1,
                   "beep" = data[, beep])

  # add observations
  df <- cbind(df, data[, y])
  # add day identifier if provided
  if (!is.null(overnight_lags)) {
    df$dayno = data[, overnight_lags]
  }

  # add between-person variables if provided
  if (!is.null(pred_random) | !is.null(outcome)) {
    out_pred_b <- out_predictors[!(out_predictors %in% c("mu", "ar", "logv"))]
    if (length(out_pred_b) == 0) {
      df <- cbind(df, data[, c(pred_random, outcome)])
      out_pred_b = NULL
    } else {
      df <- cbind(df, data[, c(pred_random, outcome, out_pred_b)])
    }
  }

  # make sure that old variable names are removed
  # in case old column names were passed to df
  if (length(y) > 1) {y_col <- paste0("y_", 1:length(y))} else {y_col <- "y"}
  colnames(df) <- c(
    "id", "dayno", "beep",
    y_col, pred_random, outcome, out_pred_b
  )

  # handling of missing values
  if (miss_handling == "remove") {
    # remove all missings
    df <- dplyr::filter(df, dplyr::if_all(dplyr::everything(), ~ !is.na(.)))
  } else if (miss.handling == "impute") {
    # set all missings to -99
    df <- dplyr::mutate(
      df,
      dplyr::across(.cols = grepl(pattern = "y", colnames(.)),
                    ~ ifelse(is.na(.), -99, .))
    )
  }
  # ensure that the id-variable is running from 1 to N
  df$id_running = as.numeric(factor(df$id))

  # ensure that the beep variable is running from 1 to TP for each subject
  df <- df %>%
    dplyr::group_by(id_running) %>%
    dplyr::mutate(beep_running = 1:dplyr::n()) %>%
    # add a running beep number per id and dayno
    dplyr::group_by(id_running, dayno) %>%
    dplyr::mutate(beep_day_running = 1:dplyr::n()) %>%
    # arrange data frame
    dplyr::arrange(id, dayno, beep_running)

  # create parameters to pass to stan
  df <- df %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      # index of each observation
      pos = 1:dplyr::n(),
      # indicate observations to include (exclude first observation per day)
      include = ifelse(beep_day_running == 1, 0, 1),
      # running id of included observations
      pos_mus = cumsum(include)
    ) %>%
    # sum of included obseervations per subject
    dplyr::group_by(id_running) %>%
    dplyr::mutate(N_obs_id = sum(include)) %>%
    dplyr::ungroup() %>%
    # arrange data frame
    dplyr::arrange(id_running, beep_running)

  # store obs per id in vector
  N_obs_id <- df %>%
    dplyr::select(id_running, N_obs_id) %>%
    dplyr::distinct() %>%
    dplyr::select(N_obs_id) %>%
    unlist()
  # maximum number of included obs for any subject
  N_obs_id_max <- max(N_obs_id)
  # maximum number of subjects
  N <- max(df$id_running)
  # maximum number of beeps per subject
  TP <- max(df$beep_running)
  # number of observations
  N_obs <- nrow(df)

  # create empty arrays to store position indices on vector of observations
  pos_mus <- array(dim = c(N, N_obs_id_max))
  pos_t <- array(dim = c(N, N_obs_id_max))
  pos_lag <- array(dim = c(N, N_obs_id_max))
  pos_mus_lag <- array(dim = c(N, N_obs_id_max))

  # fill arrays
  for (pp in 1:N) {
    pos_mus[pp, 1:N_obs_id[pp]] <-
      df$pos_mus[df$id_running == pp &df$include != 0]
    pos_t[pp, 1:N_obs_id[pp]] <-
      df$pos[df$id_running == pp & df$include != 0]
    pos_lag[pp, 1:N_obs_id[pp]] <-
      df$pos[df$id_running == pp & df$beep_running <= N_obs_id[pp]]
    pos_mus_lag[pp, 1:N_obs_id[pp]] <-
      df$pos_mus[df$id_running == pp & df$beep_running <= N_obs_id[pp]]
  }

  # long versions for latent models
  pos_mus_long <- as.vector(t(pos_mus))[as.vector(t(pos_mus))!= 0]
  pos_t_long <- as.vector(t(pos_t))[as.vector(t(pos_t))!= 0]
  pos_lag_long <- as.vector(t(pos_lag))[as.vector(t(pos_lag))!= 0]
  N_index_mus <- df$id_running[df$include != 0]
  N_index <- df$id_running

  # if external predictors are provided:
  # create matrix of person parameters to use as predictors of individual
  # parameters
  # if none are added only create a matrix of 1s for estimation
  # of intercepts (population parameter means)
  if (is.null(pred_random)) {
    n_cov <- 1
    W <- matrix(ncol = 1, nrow = N, data = 1)
  } else {
    n_cov <- 1 + length(pred_random)
    Ws <- df %>% dplyr::select(dplyr::all_of(pred_random)) %>% dplyr::distinct()
    W <- matrix(ncol = 1 + length(pred_random), nrow = N, data = 1)
    for (i in 1:length(pred_random)) {
      W[, i + 1] <- unname(unlist(Ws[, i]))
    }
  }

  # if external outcome is provided:
  # create matrix of person parameters to use as external criterion
  if (is.null(outcome)) {
    n_out <- 0
    n_out_pred <- 0
    n_out_pred_b <- 0
    out_pred_which <- matrix(nrow = 1, ncol = 0)
    out_pred_bs <- matrix(nrow = N, ncol = 0)
    out <- matrix(ncol = N, nrow = 0, data = 0)
    out_is_std <- 0
  } else {
    # outcome predictors
    out_pred_w <- out.predictors[out.predictors %in% c("mu", "ar", "logv")]
    n_out_pred <- ifelse(standardize_out == 1,
                        length(out.predictors),
                        length(out.predictors) + 1)
    out_pred_which <- sort(
      as.numeric(
        factor(out_pred_w, levels = c("mu", "ar", "logv"), labels = 1:3)
      )
    )
    out_pred_which <- matrix(
      nrow = 1, ncol = length(out_pred_w), data = out_pred_which
    )
    out_pred_b <- out.predictors[!(out.predictors %in% c("mu", "ar", "logv"))]
    n_out_pred_b <- length(out_pred_b)
    out_pred_bs <- matrix(nrow = N, ncol = n_out_pred_b)

    if (length(out_pred_b) > 0) {
      out_pred_val <- df %>%
        dplyr::select(all_of(out_pred_b)) %>%
        dplyr::distinct()
      for (i in 1:n_out_pred_b) {
        out_pred_bs[, i] <- unname(unlist(out_pred_val[, i]))
        if (standardize_out == T) {
          out_pred_bs[, i] <- scale(unlist(out_pred_val[, i]))
        }
      }
    }

    # outcome
    n_out <- length(outcome)
    outs <- df %>% dplyr::select(all_of(outcome)) %>% dplyr::distinct()
    out <- matrix(nrow = n_out, ncol = N, data = 0)
    for (i in 1:n_out) {
      out[i, ] <- unname(unlist(outs[, i]))
      if (standardize_out == T) {
        out[i, ] <- scale(unlist(outs[, i]))
      }
    }
    out_is_std <- sum(standardize_out)
  }


  # if imputation of missing values is requested:
  # fill with -99 and let Stan impute values from model
  if (miss_handling == "impute" & length(y) == 1) {
    N_miss <- sum(df$y == -99)
    pos_miss <- array(dim = c(1, N_miss), data = which(df$y == -99))
    N_miss_max <- N_miss
  } else if (miss_handling == "impute" & length(y) > 1) {
    for (ii in 1:length(y)) {
      N_miss[ii] <- sum(df[, y_col[ii]] == -99)
    }
    N_miss_max <- max(N_miss)
    pos_miss <- array(data = 0, dim = c(length(y), N_miss_max))
    for (ii in 1:length(y)) {
      pos_miss[ii, 1:N_miss[ii]] <- which(df[, y_col[ii]] == -99)
    }
  } else {
    # NEEDS UPDATING LATER
    N_miss <- array(data = 0, dim = length(y))
    N_miss_max <- array(data = 0, dim = length(y))
    pos_miss <- matrix(nrow = length(y), ncol = 0)
  }


  if(length(y) == 1){
    y_array <- df$y
  } else {
    y_array <- array(data = NA, dim = c(length(y), ncol = N_obs))
    for(ii in 1:length(y)){
      y_array[ii, ] <- unname(unlist(df[, y_col[ii]]))
    }
  }

  # create final object that can be passed to Stan
  stan_data = list(
    y = y_array,
    N = N,
    TP = TP,
    N_obs = N_obs,
    N_obs_id = N_obs_id,
    N_miss = N_miss,
    N_miss_max = N_miss_max,
    N_use = sum(df$include != 0),
    N_index_mus = N_index_mus,
    N_index = N_index,
    # N_ind = length(ts),
    N_ind = nrow(y_array)[[1]],
    n_out = n_out,
    n_out_pred = n_out_pred,
    out_pred_which = out_pred_which,
    out = out,
    out_is_std = out_is_std,
    n_out_pred_b = n_out_pred_b,
    out_pred_b = out_pred_bs,

    pos_t = pos_t,
    pos_t_long = pos_t_long,
    pos_lag = pos_lag,
    pos_lag_long = pos_lag_long,
    pos_mus = pos_mus,
    pos_mus_lag = pos_mus_lag,
    pos_mus_long = pos_mus_long,
    pos_miss = pos_miss,
    logv_is_random = ifelse(random_innovations == TRUE, 1, 0),
    n_cov = n_cov,
    W = W
    # latent case
    # N_ind = N_ind
  )

  # add Mplus data if requested
  if (add_mplus_data == TRUE) {
    mplus <- df %>%
      dplyr::select(id_running, beep_running, y_col, pred_random, outcome)
    stan_data$mplus = mplus
  }

  return(stan_data)
}

