update_model_REcors <- function(model) {
  # a helper function to call after changes were made to the number of labels
  # of the random effects in the model
  # --> update the random effect correlations included in the model

  # remove
  model <- model[model$Type != "RE correlation", ]

  # update random effect correlations
  rand.pars <- (model[model$Type == "Fixed effect" & model$isRandom == 1, "Param"])
  n_rand <- length(rand.pars)
  btw.cov_pars <- c()
  if (n_rand > 1) {
    n_cors <- (n_rand * (n_rand - 1)) / 2
    qs <- c()
    ps <- c()
    for (i in 1:(n_rand - 1)) {
      qs <- c(qs, rep(rand.pars[i], each = n_rand - i))
    }
    for (i in 2:n_rand) {
      ps <- c(ps, rep(rand.pars[i:n_rand], 1))
    }

    btw.cov_pars <- paste0("r_", qs, ".", ps)

    ## random effect correlations
    REcors <- data.frame(
      "Model" = "Structural",
      "Level" = "Between",
      "Type" = rep("RE correlation", n_cors),
      "Param" = btw.cov_pars,
      "Param_Label" = "RE Cor",
      "isRandom" = 0
    )
    model <- dplyr::bind_rows(model, mlts_model_priors(REcors, default = TRUE))
  } else if (n_rand == 1) {
    model <- model[model$Type != "RE Cor", ]
  }


  return(model)
}


replace_model_row <- function(model, row, replacement) {
  # a helper function to select a specific row in the model and replace with
  # one (or multiple) row(s) while keeping the original order

  if (row == 1) {
    model <- dplyr::bind_rows(replacement, model[2:nrow(model), ])
  } else if (row == nrow(model)) {
    model <- dplyr::bind_rows(model[1:(nrow(model) - 1), ], replacement)
  } else {
    model <- dplyr::bind_rows(
      model[1:(row - 1), ],
      replacement,
      model[(row + 1):nrow(model), ]
    )
  }

  return(model)
}

extract_indicator_info <- function(model, level = "Within", type = "Loading", incl.pos_p = FALSE) {
  # a helper function to extract indicator information

  # create a table where all indocators are present
  # get subset
  info <- model[model$Level == level & model$Type == type, ]

  # extract infos from parameter
  inds <- unlist(lapply(info$Param, function(x) {
    strsplit(x, split = "_")[[1]][2]
  }))
  param <- unique(unlist(lapply(info$Param, function(x) {
    strsplit(x, split = "_")[[1]][1]
  })))

  ind.info <- data.frame(
    "q" = unlist(lapply(inds, function(x) {
      strsplit(x, split = ".", fixed = TRUE)[[1]][1]
    })),
    "p" = unlist(lapply(inds, function(x) {
      strsplit(x, split = ".", fixed = TRUE)[[1]][2]
    }))
  )

  # add general indicator number
  if (incl.pos_p == TRUE) {
    ind.info$p_pos <- 1:nrow(ind.info)
  }
  # add parameter
  if(level== "Within" & type == "Loading"){
    ind.info[, paste0(param, "_isFree")] <- ifelse(info$Constraint != "= 1" , 1, 0)
  } else {
    ind.info[, paste0(param, "_isFree")] <- ifelse(info$Constraint == "free", 1, 0)
  }
  # add loading parameter constraints
  if(type == "Loading"){
    ind.info[, paste0(param, "_isEqual")] <- info$Constraint
  }

  return(ind.info)
}

# change colnames in summary function
change_colnames <-  function(data, cols) {
  names <- colnames(data)
  names[grepl("Param", colnames(data))] <- ""
  names[grepl("50%", colnames(data))] <- "Median"
  names[grepl("mean", colnames(data))] <- "Mean"
  names[grepl("sd", colnames(data))] <- "SD"
  return(names)
}

# function to evaluate input of t0_effects
eval_t0_effects <- function(t0_input, q){

  # initial check
  check1 <- sum(startsWith(prefix = "phi(0)_", x = t0_input))
  if(check1 != length(t0_input)){
    stop("Invalid input of 'incl_t0_effects', see ?mlts_model.")
  }

  t0_effs <- lapply(t0_input, function(x){
    df <- data.frame(
      "DV" = strsplit(strsplit(x, split = "_")[[1]][2], split = "")[[1]][1],
      "IV" = strsplit(strsplit(x, split = "_")[[1]][2], split = "")[[1]][2]
    )
    df
  })
  t0_effs <- do.call(rbind, t0_effs)

  # post sanity checks
  # bidirectional effects
  n_bi <- paste0(t0_effs$DV,t0_effs$IV)  %in% paste0(t0_effs$IV,t0_effs$DV)

  if(sum(n_bi)>0){
    stop("Invalid input of 'incl_t0_effects': bidrectional paths are not allowed (e.g., phi(0)_21 and phi(0)_12.")
  }

  if(max(t0_effs$DV) > q | max(t0_effs$IV) > q){
    stop("Invalid input of 'incl_t0_effects': input refers to variables outside the number of included constructs ('q').")
  }


  return(t0_effs)
}


# function to evaluate input of incl_interaction_effects
eval_int_effects <- function(int_input, q){

  # initial input check
  checks = list()
  checks[[1]] <- sum(startsWith(prefix = "phi(i)_", x = int_input))
  checks[[2]] <- sum(unlist(lapply(int_input, function(x){substr(x, 9,9) == "."})))
  checks[[3]] <- sum(unlist(lapply(int_input, function(x){substr(x, 11,11) == "("})))
  checks[[4]] <- sum(unlist(lapply(int_input, function(x){substr(x, 13,13) == ")"})))
  checks[[5]] <- sum(unlist(lapply(int_input, function(x){substr(x, 15,15) == "("})))
  checks[[6]] <- sum(unlist(lapply(int_input, function(x){substr(x, 17,17) == ")"})))
  if(any(unlist(checks) != length(int_input))){
    stop("Invalid input of 'incl_interaction_effects', see ?mlts_model.")
  }

  int_effs <- lapply(int_input, function(x){
    # remove prefix
    df <- data.frame(
      "Param" = x,
      "DV"  = substr(x, 8,8),
      "IV1" = substr(x, 10,10),
      "IV1lag" = substr(x, 12,12),
      "IV2" = substr(x, 14,14),
      "IV2lag" = substr(x, 16,16)
    )
    df
  })
  int_effs <- do.call(rbind, int_effs)

  # post sanity checks
  # main effects included?

  # all variables in model?
  if(max(as.integer(int_effs$DV)) > q |
     max(as.integer(int_effs$IV1)) > q |
     max(as.integer(int_effs$IV2)) > q){
    stop("Invalid input of 'incl_t0_effects': input refers to variables outside the number of included constructs ('q').")
  }


  return(int_effs)
}
