#' Correlations for time series variables
#'
#' Computes correlations between time series variables at the between-person
#' or within-person level.
#'
#' At the between-person level, correlations are computed based on person means.
#' At the within-person level, correlations are computed within clusters and
#' can optionally be aggregated across clusters using a summary function.
#'
#' @param data A data.frame containing the time series data.
#' @param ts Character vector of column names representing the time series variables.
#'   Must contain at least two variables.
#' @param id Character string specifying the subject identifier column.
#' @param type Character string specifying the level of analysis.
#'   Must be either \code{"between"} or \code{"within"}.
#' @param as_mat Logical. If \code{TRUE}, returns correlation matrices.
#'   If \code{FALSE}, returns a data.frame with pairwise correlations and confidence intervals.
#' @param within_agg_fn Optional function used to aggregate within-person correlations
#'   across individuals (e.g., \code{mean}, \code{median}). Only used when
#'   \code{type = "within"} and \code{as_mat = TRUE}.
#'
#' @return
#' If \code{type = "between"}:
#' \itemize{
#'   \item If \code{as_mat = TRUE}: A correlation matrix.
#'   \item If \code{as_mat = FALSE}: A data.frame with variables, correlations,
#'   and confidence intervals.
#' }
#'
#' If \code{type = "within"}:
#' \itemize{
#'   \item If \code{as_mat = TRUE} and \code{within_agg_fn = NULL}: A list of
#'   individual correlation matrices.
#'   \item If \code{as_mat = TRUE} and \code{within_agg_fn} is provided:
#'   A single aggregated correlation matrix.
#'   \item If \code{as_mat = FALSE}: A data.frame with individual-level
#'   correlations and confidence intervals.
#' }
#'
#' @details
#' Correlations are computed using Pearson's correlation coefficient with
#' pairwise complete observations. Confidence intervals are obtained using
#' \code{\link[stats]{cor.test}}. If estimation fails, \code{NA} values are returned.
#'
#' @examples
#' \dontrun{
#' # Specify a three-indicator AR(1) model
#' model <- mlts_model(q = 1, p = 3)
#'
#' # Generate data under the specified model
#' df <- mlts_sim(model = model, N = 50, TP = 100, default = T)$data
#'
#' # Between-person correlations
#' mlts_cor_ts(data = df, ts = c("Y1.1", "Y1.2", "Y1.3"), id = "ID")
#'
#' # Within-person correlations (list of matrices)
#' mlts_cor_ts(df, ts = c("Y1.1", "Y1.2", "Y1.3"), id = "ID", type = "within")
#'
#' # Aggregated within-person correlations
#' mlts_cor_ts(df, ts = c("Y1.1", "Y1.2", "Y1.3"), id = "ID",
#'             type = "within", within_agg_fn = mean)
#' }
#'
#' @export

mlts_cor_ts <- function(data, ts, id, type = "between", as_mat = TRUE, within_agg_fn = NULL){
  type <- match.arg(type, choices = c("between", "within"))
  if (length(ts) == 1){
    stop("ts must contain at least two items to calculate a correlation", call. = FALSE)
  }
  n_items <- length(ts)

  # between level correlations
  if(type == "between"){
    person_means <- stats::aggregate(x = data[, ts ], by = list(data[[id]]), FUN = mean, na.rm = TRUE)
    cor_mat <- stats::cor(person_means[,-1], method = "pearson", use = "pairwise.complete.obs")

    if (as_mat){
      return(cor_mat)

    } else {

      ts_vars_df = data.frame(t(utils::combn(ts, 2)))

      results = lapply(1:nrow(ts_vars_df), function(x){
        item = ts_vars_df$X1[x]
        item2 = ts_vars_df$X2[x]

        c_test <- tryCatch({
          stats::cor.test(person_means[, item],
                          person_means[, item2], method = "pearson")
        }, error = function(e) {
          list(conf.int = c(NA, NA))
        })

        temp_df <- data.frame(
          Var1 = ts[item],
          Var2 = ts[item2],
          r = cor_mat[item, item2],
          CI.lb = c_test$conf.int[1],
          CI.ub = c_test$conf.int[2]
        )

      })

      return(do.call(rbind, results))

    }
  }


  if (type == "within"){

    u_ids <- unique(data[[id]])
    matrices <- list()
    for ( i in u_ids){
      person_data <- data[data[[id]] == i, ts]
      matrices[[as.character(i)]] <- stats::cor(person_data, method = "pearson", use = "pairwise.complete.obs")
    }

    if(as_mat){

      if (!is.null(within_agg_fn)){
        temp_cors <- numeric(length(matrices))
        final_mat <- matrix(nrow = n_items, ncol = n_items)
        colnames(final_mat) <- ts
        rownames(final_mat) <- ts

        for (item in 1:(n_items-1)){
          for (item2 in (item+1):n_items){
            for ( mat in seq_along(matrices)){
              temp <- matrices[[mat]]
              temp_cors[mat] <- temp[item, item2]
            }
            fn_value <- within_agg_fn(temp_cors, na.rm = TRUE)
            final_mat[item, item2] <- fn_value
            final_mat[item2, item] <- fn_value
          }
        }
        diag(final_mat) <- 1
        return(final_mat)
      }
      return(matrices)
    }
    else {
      results <- list()
      list_counter = 1
      for (i in u_ids){
        person_data <- data[data[[id]] == i, ts]
        for (item in 1:(n_items-1)){
          for (item2 in (item+1):n_items){
            vec1 <- person_data[, item]
            vec2 <- person_data[, item2]

            c_test <- tryCatch({
              stats::cor.test(vec1,vec2, method = "pearson")
            }, error = function(e){
              list(estimate = NA, conf.int = c(NA, NA))
            })

            temp_df <- data.frame(
              ID = i,
              Var1 = ts[item],
              Var2 = ts[item2],
              r = c_test$estimate,
              CI.lb = c_test$conf.int[1],
              CI.ub = c_test$conf.int[2]
            )
            results[[list_counter]] <- temp_df
            list_counter <- list_counter +1
          }
        }
      }
      out = do.call(rbind, results)
      rownames(out) <- NULL
      return(out)
    }
  }

}
