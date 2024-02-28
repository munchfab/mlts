update_model_REcors <- function(VARmodel) {
  # a helper function to call after changes were made to the number of labels
  # of the random effects in the model
  # --> update the random effect correlations included in the VARmodel

  # remove
  VARmodel <- VARmodel[VARmodel$Type != "RE correlation", ]

  # update random effect correlations
  rand.pars <- (VARmodel[VARmodel$Type == "Fix effect" & VARmodel$isRandom == 1, "Param"])
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
    VARmodel <- plyr::rbind.fill(VARmodel, mlts_model_priors(REcors, default = T))
    # consider dplyr because plyr is deprecated
    # VARmodel <- dplyr::bind_rows(VARmodel, mlts_model_priors(REcors, default = T))
  } else if (n_rand == 1) {
    VARmodel <- VARmodel[VARmodel$Type != "RE Cor", ]
  }


  return(VARmodel)
}


replace_model_row <- function(VARmodel, row, replacement) {
  # a helper function to select a specific row in the model and replace with
  # one (or multiple) row(s) while keeping the original order

  if (row == 1) {
    VARmodel <- plyr::rbind.fill(replacement, VARmodel[2:nrow(VARmodel), ])
    # consider dplyr because plyr is deprecated
    # VARmodel <- dplyr::bind_rows(replacement, VARmodel[2:nrow(VARmodel), ])
  } else if (row == nrow(VARmodel)) {
    VARmodel <- plyr::rbind.fill(VARmodel[1:(nrow(VARmodel) - 1), ], replacement)
    # consider dplyr because plyr is deprecated
    # VARmodel <- dplyr::bind_rows(VARmodel[1:(nrow(VARmodel) - 1), ], replacement)
  } else {
    VARmodel <- plyr::rbind.fill(
      VARmodel[1:(row - 1), ],
      replacement,
      VARmodel[(row + 1):nrow(VARmodel), ]
    )
    # consider dplyr because plyr is deprecated
    # VARmodel <- dplyr::bind_rows(
    #   VARmodel[1:(row - 1), ],
    #   replacement,
    #   VARmodel[(row + 1):nrow(VARmodel), ]
    # )
  }

  return(VARmodel)
}

extract_indicator_info <- function(VARmodel, level = "Within", type = "Loading", incl.pos_p = F) {
  # a helper function to extract indicator information

  # create a table where all indocators are present
  # get subset
  info <- VARmodel[VARmodel$Level == level & VARmodel$Type == type, ]

  # extract infos from parameter
  inds <- unlist(lapply(info$Param, function(x) {
    strsplit(x, split = "_")[[1]][2]
  }))
  param <- unique(unlist(lapply(info$Param, function(x) {
    strsplit(x, split = "_")[[1]][1]
  })))

  ind.info <- data.frame(
    "q" = unlist(lapply(inds, function(x) {
      strsplit(x, split = ".", fixed = T)[[1]][1]
    })),
    "p" = unlist(lapply(inds, function(x) {
      strsplit(x, split = ".", fixed = T)[[1]][2]
    }))
  )

  # add general indicator number
  if (incl.pos_p == T) {
    ind.info$p_pos <- 1:nrow(ind.info)
  }
  # add parameter
  ind.info[, paste0(param, "_isFree")] <- ifelse(info$Constraint == "free", 1, 0)

  return(ind.info)
}
