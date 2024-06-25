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
  ind.info[, paste0(param, "_isFree")] <- ifelse(info$Constraint == "free", 1, 0)

  return(ind.info)
}

# change colnames in summary function
change_colnames <-  function(data, cols) {
  names <- colnames(data)
  names[grepl("Param", colnames(data))] <- ""
  names[grepl("50%", colnames(data))] <- "Post. Median"
  names[grepl("mean", colnames(data))] <- "Post. Mean"
  names[grepl("sd", colnames(data))] <- "Post. SD"
  return(names)
}
