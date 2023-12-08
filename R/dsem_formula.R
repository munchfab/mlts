#' Title
#'
#' @param y The variable in `data` used for the autoregressive process
#' (should be the same as the response variable).
#' @param p A non-negative integer specifying the autoregressive (AR) order.
#' Default is 1.
#' @inheritParams dsem
#' @return
#' @export
#'
#' @examples
ar <- function(y, p = 1) {
  ar_call <- paste0("ar(y, p = ", p, ")")
  return()
}


#' Title
#'
#' @param formula
#'
#' @return
#' @export
#'
#' @examples
dsem_formula <- function(formula) {

  # old syntax ================================================================
  # dsem_terms <- terms(as.formula(formula))
  # mf <- model.frame(dsem_terms, data = data)
  # y <- model.response(mf, type = "numeric")
  # all_terms <- attr(dsem_terms, "term.labels")

  # extract right hand side of formula
  # rhs <- rlang::f_rhs(formula)
  # all_terms <- deparse(rhs)
  # extract all terms with brackets (random effect formulas and special syntax)
  # bracket_terms <- regmatches(
  #   all_terms,
  #   gregexpr("\\w*(?=\\().+?(?<=\\))", all_terms, perl = TRUE)
  # )
  # extract random effect terms
  # re_terms <- gsub(
  #   pattern = "\\s", # delete all space characters
  #   replacement = "",
  #   bracket_terms[[1]][grepl("\\|", bracket_terms[[1]])]
  # )
  # specials <- gsub(
  #   pattern = "\\s", # delete all space characters
  #   replacement = "",
  #   bracket_terms[[1]][!grepl("\\|", bracket_terms[[1]])]
  # )

  # new syntax ================================================================
  # extract response from formula
  response <- all.vars(formula)[1]
  # extract all terms from formula
  all_terms <- labels(terms(formula))
  # extract fixed effect terms
  fe_terms <- gsub(
    pattern = "\\s", # delete all space characters
    replacement = "",
    all_terms[!grepl("\\|", all_terms)]
  )
  # extract random effect terms from all_terms
  re_terms <- gsub(
    pattern = "\\s", # delete all space characters
    replacement = "",
    all_terms[grepl("\\|", all_terms)]
  )
  # extract special syntax
  specials <- gsub(
    pattern = "\\s", # delete all space characters
    replacement = "",
    all_terms[grepl("\\w*\\(", all_terms)]
  )
  # extract components for Stan model
  beep <- regmatches(
    re_terms, gregexpr("(?<=\\|)\\w*(?=\\/)", re_terms, perl = TRUE)
  )[[1]]
  # beep <- regmatches(
  #   re_terms, gregexpr("(?<=\\|).+?(?=\\/)", re_terms, perl = TRUE)
  # )[[1]]
  id <- regmatches(
    re_terms, gregexpr("(?<=\\/)\\w*", re_terms, perl = TRUE)
  )[[1]]

  # store in list
  dsem_terms <- rstan::nlist(
    all_terms, fe_terms, re_terms, specials,
    beep, id, response
  )

  return(dsem_terms)
}
