#' Create Path Diagrams from VARmodel object
#'
#' @param VARmodel The VARmodel object.
#' @param data An optional `data.frame` including the variables to be used in
#' the path diagram.
#' @param labels An optional character string inluding the names to be used in
#' the path diagram.
#' @param add.png Logical. Set to `TRUE` to transform created PDF as .png files
#'  using `pdftools::pdf_convert`.
#' @return An Rmarkdown file that is automatically rendered to a pdf document.
#' @export
#'
#' @examples
VARmodelPaths <- function(VARmodel, data = NULL, labels = NULL, add.png = FALSE) {


  # extract model infos
  infos <- VARmodelEval(VARmodel)

  # extract model data frame
  model <- VARmodel

  # create empty markdown file, delete if already existing
  if (file.exists("pathmodel.rmd")) {
    file.remove("pathmodel.rmd")
  }
  rmarkdown::draft(file = "pathmodel.rmd",
                   template = "pathmodel",
                   package = "dsemr",
                   edit = FALSE)

  # latex figure begin and end
  begin_figure <- "\n\\begin{figure}\n\\centering"
  end_figure <- "\\end{figure}\n"

  # create empty string with tikz specifications
  begin_tikz <- "\\begin{tikzpicture}[
    auto, > = latex, align=center,
  	latent/.style = {
  	  circle, draw, thick, inner sep = 2pt, minimum width = \\Radius
  	},
  	manifest/.style = {
  	  rectangle, draw, thick, inner sep = 0pt, minimum size = \\Radius/2
  	},
  	intercept/.style = {
  	  regular polygon,regular polygon sides = 3, draw, thick,
  	  inner sep = 0pt, minimum size = \\Radius
  	},
  	mean/.style = {
  	  regular polygon, regular polygon sides = 3, draw, thick,
  	  inner sep = 0pt, minimum size = 8mm
  	},
  	path/.style = {
  	  arrows = ->, thick, > = {stealth[]}
  	},
  	error/.style = {
  	  circle, draw = none, fill = none, thick,
  	  inner sep = 0pt, minimum size = 5mm
  	},
  	var/.style = {
  	  <->, thick, > = {stealth[]}, bend right = 270, looseness = 2
  	},
  	cov/.style = {
  	  <->, thick, > = {stealth[]}, bend right = 300
  	},
  	% style to add a circle in the middle of a path
    random/.style = {
      postaction = {
        decorate, decoration = {
          markings,
          mark = at position .5 with {
            \\draw[fill = black] circle[radius = 2pt];
          }
        }
      }
    },
    random_cl/.style = {
      postaction = {decorate, decoration = {
        markings,
        mark = at position .2 with {
          \\draw[fill = black] circle[radius = 2pt];}
        }
      }
    },
  	]"

  # end tikz picture
  end_tikz <- "\\end{tikzpicture}\n"


  # decomposition #############################################################

  dc_caption <- paste0(
    "% caption\n\\node  [above = 1em, align = center]  at  ",
    "(current bounding box.north)  {Decomposition.};"
  )

  # one manifest construct
  if (infos$q == 1) {
    dc <- "
    % draw decomposition
    \\node  [manifest]  (y1t)  {$y_{1,t}$};
    \\node  [latent]  (y1wt)  [above = 2.5em of y1t]  {$y_{1,t}^w$};
    \\node  [latent]  (mu1)  [below = 2.5em of y1t]  {$\\mu_{1}$};

    % draw paths
    \\draw  [path]  (y1wt)  to node  []  {}  (y1t);
    \\draw  [path]  (mu1)  to node  []  {}  (y1t);
    "
    if (infos$p > 1) {
      # initiate empty vectors to store tikz nodes
      ind_vec <- c() # indicator variables
      # latent variables (within and between)
      wlat <- "\\node  [latent]  (eta1wt)  [above = 4em of c]  {$\\eta^w_{1,t}$};"
      blat <- "\\node  [latent]  (eta1b)  [below = 4em of c]  {$\\eta^b_{1}$};"
      wlat_paths_vec <- c() # paths from latent variables to indicators
      blat_paths_vec <- c() # paths from latent variables to indicators
      epswt_vec <- c() # residuals within
      epswt_paths_vec <- c() # paths from residuals within to indicators
      epsb_vec <- c() # residuals within
      epsb_paths_vec <- c() # paths from residuals within to indicators
      for (i in 1:infos$p) {
        if (i == 1) {
          # first indicator variable
          ind_vec <- c(ind_vec, paste0(
            "\\node  [manifest]  (y", i, "t)  {$y_{", i, ",t}$};"
          ))
          # first error (labeled)
          epswt_vec <- c(epswt_vec, paste0(
            "\\node  [error]  (eps", i, "wt)  [above left = 1em of y",
            i, "t]  {\\scriptsize$\\varepsilon^w_{", i, ",t}$};"
          ))
          epsb_vec <- c(epsb_vec, paste0(
            "\\node  [error]  (eps", i, "b)  [below left = 1em of y",
            i, "t]  {\\scriptsize$\\varepsilon^b_{", i, "}$};"
          ))
        } else {
          # all other indicators
          ind_vec <- c(ind_vec, paste0(
            "\\node  [manifest]  (y", i, "t)  [right = 1.5em of y",
            i - 1, "t]  {$y_{", i, ",t}$};"
          ))
          # residuals
          epswt_vec <- c(epswt_vec, paste0(
            "\\node  [error]  (eps", i, "wt)  [above left = .75em of y", i, "t]  {};"
          ))
          epsb_vec <- c(epsb_vec, paste0(
            "\\node  [error]  (eps", i, "b)  [below left = .75em of y", i, "t]  {};"
          ))
          if (i == infos$p) {
            # place coordinates between first and last indicator variables
            # of construct
            coord <- paste0(
              "\\node  [coordinate]  (c)  at ($(y1t) !0.5! (y", i, "t)$)  {};"
            )
          }
        }
        # draw paths
        wlat_paths_vec <- c(wlat_paths_vec, paste0(
          "\\draw  [path]  (eta1wt)  to node  [fill = white, anchor = center]
            {\\scriptsize", ifelse(
              # label path with constraint if necessary
              model[model$Param == paste0("lambdaW_1.", i), "Constraint"] == "= 1",
              "$1$}", paste0("$\\lambda^w_{", i, "}$}")
            ), "  (y", i, "t.north);"
        ))
        blat_paths_vec <- c(blat_paths_vec, paste0(
          "\\draw  [path]  (eta1b)  to node  [fill = white, anchor = center]
            {\\scriptsize", ifelse(
              # label path with constraint if necessary
              model[model$Param == paste0("lambdaB_1.", i), "Constraint"] == "= 1",
              "$1$}", paste0("$\\lambda^b_{", i, "}$}")
            ), "  (y", i, "t.south);"
        ))
        epswt_paths_vec <- c(epswt_paths_vec, paste0(
          "\\draw  [path]  (eps", i, "wt)  to node  []  {}  (y", i, "t.north west);"
        ))
        epsb_paths_vec <- c(epsb_paths_vec, paste0(
          "\\draw  [path]  (eps", i, "b)  to node  []  {}  (y", i, "t.south west);"
        ))
      }
      # paste together
      ind <- paste(ind_vec, collapse = "\n")
      wlat_paths <- paste(wlat_paths_vec, collapse = "\n")
      blat_paths <- paste(blat_paths_vec, collapse = "\n")
      epswt <- paste(epswt_vec, collapse = "\n")
      epswt_paths <- paste(epswt_paths_vec, collapse = "\n")
      epsb <- paste(epsb_vec, collapse = "\n")
      epsb_paths <- paste(epsb_paths_vec, collapse = "\n")
      dc <- paste(ind, coord, wlat, blat,
                  wlat_paths, blat_paths,
                  epswt, epswt_paths, epsb, epsb_paths, sep = "\n")
    }
  } else { # for q > 1
    dc <- paste0(
    "
    % draw decomposition
    \\node  [manifest]  (y1t)  {$y_{1,t}$};
    \\node  [latent]  (y1wt)  [above = 2.5em of y1t]  {$y_{1,t}^w$};
    \\node  [latent]  (mu1)  [below = 2.5em of y1t]  {$\\mu_{1}$};

    % draw nodes
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    infos$q, "}
    \\node  [manifest]  (y\\i t)  [right = 2.5em of y\\lasti t]  {$y_{\\i,t}$};

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    infos$q, "}
    \\node  [latent]  (y\\i wt)  [above = 2.5em of y\\i t]  {$y_{\\i,t}^w$};

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ...,",
    infos$q, "}
    \\node  [latent]  (mu\\i)  [below = 2.5em of y\\i t]  {$\\mu_{\\i}$};

    % draw paths
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
    infos$q, "}
    \\draw  [path]  (y\\i wt)  to node  []  {}  (y\\i t);

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
    infos$q, "}
    \\draw  [path]  (mu\\i)  to node  []  {}  (y\\i t);"
    )
    if (any(infos$p > 1)) {
      # initiate empty vectors to store tikz nodes
      ind_vec <- c() # indicator variables
      coord_vec <- c() # coordinates for positioning of latent variables
      wlat_vec <- c() # within-level latent variables
      blat_vec <- c() # between-level latent variables
      wlat_paths_vec <- c() # paths from latent variables to indicators
      blat_paths_vec <- c() # paths from latent variables to indicators
      epswt_vec <- c() # residuals within
      epswt_paths_vec <- c() # paths from residuals within to indicators
      epsb_vec <- c() # residuals within
      epsb_paths_vec <- c() # paths from residuals within to indicators
      for (i in 1:infos$q) {
        wlat_vec <- c(wlat_vec, paste0(
          "\\node  [latent]  (etaw", i, "t)  [above = 4em of c", i,
          "]  {$\\eta^w_{", i, ",t}$};"
        ))
        blat_vec <- c(blat_vec, paste0(
          "\\node  [latent]  (etab", i, ")  [below = 4em of c", i,
          "]  {$\\eta^b_{", i, "}$};"
        ))
        for (j in 1:infos$p[i]) {
          if (i == 1 & j == 1) {
            # first indicator variable
            ind_vec <- c(ind_vec, paste0(
              "\\node  [manifest]  (y", i, j, "t)  {$y_{", i, j, ",t}$};"
            ))
            # first error (labeled)
            epswt_vec <- c(epswt_vec, paste0(
              "\\node  [error]  (eps", i, j, "wt)  [above left = 1em of y",
              i, j, "t]  {\\scriptsize$\\varepsilon^w_{", i, j, ",t}$};"
            ))
            epsb_vec <- c(epsb_vec, paste0(
              "\\node  [error]  (eps", i, j, "b)  [below left = 1em of y",
              i, j, "t]  {\\scriptsize$\\varepsilon^b_{", i, j, "}$};"
            ))
          } else if (i > 1 & j == 1) {
            # first indicator variable of next construct
            ind_vec <- c(ind_vec, paste0(
              "\\node  [manifest]  (y", i, j, "t)  [right = 1.5em of y",
              i - 1, infos$p[i - 1], "t]  {$y_{", i, j, ",t}$};"
            ))
            epswt_vec <- c(epswt_vec, paste0(
              "\\node  [error]  (eps", i, j, "wt)  [above left = .75em of y", i, j, "t]  {};"
            ))
            epsb_vec <- c(epsb_vec, paste0(
              "\\node  [error]  (eps", i, j, "b)  [below left = .75em of y", i, j, "t]  {};"
            ))
          } else {
            # all other indicators
            ind_vec <- c(ind_vec, paste0(
              "\\node  [manifest]  (y", i, j, "t)  [right = 1.5em of y",
              i, j - 1, "t]  {$y_{", i, j, ",t}$};"
            ))
            epswt_vec <- c(epswt_vec, paste0(
              "\\node  [error]  (eps", i, j, "wt)  [above left = .75em of y", i, j, "t]  {};"
            ))
            epsb_vec <- c(epsb_vec, paste0(
              "\\node  [error]  (eps", i, j, "b)  [below left = .75em of y", i, j, "t]  {};"
            ))
            if (j == infos$p[i]) {
              # place coordinates between first and last indicator variables
              # of construct i
              coord_vec <- c(coord_vec, paste0(
                "\\node  [coordinate]  (c", i,
                ")  at ($(y", i, "1t) !0.5! (y", i, j, "t)$)  {};"
              ))
            }
          }
          # draw paths
          wlat_paths_vec <- c(wlat_paths_vec, paste0(
            "\\draw  [path]  (etaw",
            i, "t)  to node  [fill = white, anchor = center]  {\\scriptsize",
            ifelse(
              # label path with constraint if necessary
              model[model$Param == paste0("lambdaW_", i, ".", j), "Constraint"] == "= 1",
              "$1$}", paste0("$\\lambda^w_{", i, j, "}$}")
            ), "  (y", i, j, "t.north);"
          ))
          blat_paths_vec <- c(blat_paths_vec, paste0(
            "\\draw  [path]  (etab",
            i, ")  to node  [fill = white, anchor = center]  {\\scriptsize",
            ifelse(
              # label path with constraint if necessary
              model[model$Param == paste0("lambdaB_", i, ".", j), "Constraint"] == "= 1",
              "$1$}", paste0("$\\lambda^b_{", i, j, "}$}")
            ), "  (y", i, j, "t.south);"
          ))
          epswt_paths_vec <- c(epswt_paths_vec, paste0(
            "\\draw  [path]  (eps", i, j, "wt)  to node  []  {}  (y", i, j, "t.north west);"
          ))
          epsb_paths_vec <- c(epsb_paths_vec, paste0(
            "\\draw  [path]  (eps", i, j, "b)  to node  []  {}  (y", i, j, "t.south west);"
          ))
        }
      }
      # paste together
      ind <- paste(ind_vec, collapse = "\n")
      coord <- paste(coord_vec, collapse = "\n")
      wlat <- paste(wlat_vec, collapse = "\n")
      blat <- paste(blat_vec, collapse = "\n")
      wlat_paths <- paste(wlat_paths_vec, collapse = "\n")
      blat_paths <- paste(blat_paths_vec, collapse = "\n")
      epswt <- paste(epswt_vec, collapse = "\n")
      epswt_paths <- paste(epswt_paths_vec, collapse = "\n")
      epsb <- paste(epsb_vec, collapse = "\n")
      epsb_paths <- paste(epsb_paths_vec, collapse = "\n")
      dc <- paste(ind, coord, wlat, blat,
                  wlat_paths, blat_paths,
                  epswt, epswt_paths, epsb, epsb_paths, sep = "\n")
    }
  }

  # paste together
  decomposition <- paste(
    # begin_figure, # not needed
    # dc_caption,
    begin_tikz,
    dc,
    dc_caption,
    end_tikz,
    # end_figure,
    sep = "\n"
  )

  # paste decomposition to markdown
  cat(decomposition, file = "pathmodel.rmd", append = TRUE)


  # within-model ##############################################################

  # within-model caption
  wm_caption <- paste0(
    "% caption\n\\node  [above = 1em, align = center]  at  ",
    "(current bounding box.north)  {Within-model.};"
  )

  # draw nodes for q time-series constructs
  if (infos$q == 1) {
    wm <- paste0(
      "% draw within-level structural model
      \\node  [latent]  (y1wt-1)  {$y_{1,t-1}^w$};
      \\node  [latent]  (y1wt)  [right = 5em of y1wt-1]  {$y_{1,t}^w$};
      \\node  [latent]  (delta1t)  [right = 2.5em of y1wt]  {$\\delta_{1,t}$};

      \\foreach \\j [remember = \\j as \\lastj (initially 1)] in {2, ..., ",
      infos$maxLag, "}
      \\node  [latent]  (y1wt-\\j)  ",
      "[left = 5em of y1wt-\\lastj]  {$y_{1,t-\\j}^w$};

      % draw paths
      \\draw  [path]  (delta1t)  to node  []  {}  (y1wt);
      "
    )
    if (infos$p > 1) {
      wm <- paste0(
        "% draw within-level structural model
        \\node  [latent]  (eta1wt-1)  {$\\eta_{1,t-1}^w$};
        \\node  [latent]  (eta1wt)  [right = 5em of eta1wt-1]  {$\\eta_{1,t}^w$};
        \\node  [latent]  (delta1t)  [right = 2.5em of eta1wt]  {$\\delta_{{\\eta^w_{1}},t}$};

        \\foreach \\j [remember = \\j as \\lastj (initially 1)] in {2, ..., ",
        infos$maxLag, "}
        \\node  [latent]  (eta1wt-\\j)  ",
        "[left = 5em of eta1wt-\\lastj]  {$\\eta_{1,t-\\j}^w$};

        % draw paths
        \\draw  [path]  (delta1t)  to node  []  {}  (eta1wt);
      "
      )
    }
  } else { # for q > 1
    wm <- paste0(
    "
    % draw within-level structural model
    \\node  [latent]  (y1wt-1)  {$y_{1,t-1}^w$};
    \\node  [latent]  (y1wt)  [right =5em of y1wt-1]  {$y_{1,t}^w$};
    \\node  [latent]  (delta1t)  [right =2.5em of y1wt]  {$\\delta_{1,t}$};

    % draw nodes
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
      infos$q, "}
    \\node  [latent]  (y\\i wt-1)  [below = 5em of y\\lasti wt-1]  {$y_{\\i ,t-1}^w$};",

    ifelse(
      infos$maxLag != 1,
      paste0(
        "\\foreach \\j [remember = \\j as \\lastj (initially 1)] in {2, ..., ",
          infos$maxLag, "}
        \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ..., ",
          infos$q, "}
        \\node  [latent]  (y\\i wt-\\j)  [left = 5em of y\\i wt-\\lastj]  {$y_{\\i ,t-\\j}^w$};"
      ),
      paste0("")
    ),

    "\\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    infos$q, "}
    \\node  [latent]  (y\\i wt)  [below = 5em of y\\lasti wt]  {$y_{\\i ,t}^w$};

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    infos$q, "}
    \\node  [latent]  (delta\\i t)  [below = 5em of delta\\lasti t]  {$\\delta_{\\i ,t}$};

    % draw paths from residuals to yiwt
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
    infos$q, "}
    \\draw  [path]  (delta\\i t)  to node  []  {}  (y\\i wt);
    ")
    if (any(infos$p > 1)) {
      wm <- paste0(
      "
      % draw within-level structural model
      \\node  [latent]  (eta1wt-1)  {$\\eta_{1,t-1}^w$};
      \\node  [latent]  (eta1wt)  [right = 5em of eta1wt-1]  {$\\eta_{1,t}^w$};
      \\node  [latent]  (delta1t)  [right = 2.5em of eta1wt]  {$\\delta_{1,t}$};

      % draw nodes
      \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
        infos$q, "}
      \\node  [latent]  (eta\\i wt-1)  [below = 5em of eta\\lasti wt-1]  {$\\eta_{\\i ,t-1}^w$};",

      ifelse(
        infos$maxLag != 1,
        paste0(
          "\\foreach \\j [remember = \\j as \\lastj (initially 1)] in {2, ..., ",
            infos$maxLag, "}
          \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ..., ",
            infos$q, "}
          \\node  [latent]  (eta\\i wt-\\j)  [left = 5em of eta\\i wt-\\lastj]  {$\\eta_{\\i ,t-\\j}^w$};"
        ),
        paste0("")
      ),

      "\\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
        infos$q, "}
      \\node  [latent]  (eta\\i wt)  [below = 5em of eta\\lasti wt]  {$\\eta_{\\i ,t}^w$};

      \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
        infos$q, "}
      \\node  [latent]  (delta\\i t)  [below = 5em of delta\\lasti t]  {$\\delta_{\\i ,t}$};

      % draw paths from residuals to etaiwt
      \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
        infos$q, "}
      \\draw  [path]  (delta\\i t)  to node  []  {}  (eta\\i wt);
      ")
    }
  }

  # store number of phi-parameters for loops
  all_phis <- model[grepl("phi", model$Param) & grepl("Fix", model$Type), ]
  all_phis$names <- gsub(
    # replace underscore digit with latex subscript
    "(\\w+)(\\(\\d\\))_(\\d+)", "$\\\\\\1_{\\3;\\2}$", all_phis$Param
  )
  # generate path start and end
  all_phis$from <- gsub(
    "(\\w+)\\((\\d)\\)_(\\d)(\\d)", "\\3wt-\\2", all_phis$Param
  )
  all_phis$to <- gsub(
    "(\\w+)\\((\\d)\\)_(\\d)(\\d)", "\\4wt", all_phis$Param
  )
  # repeat AR(1) effects if maxLag > 1
  if (infos$maxLag > 1) {
    rep_phis_list <- list()
    for (i in 1:(infos$maxLag - 1)) {
      rep_phis_list[[i]] <- all_phis[grep("phi\\(1\\)", all_phis$Param), ]
      rep_phis_list[[i]]$from <- gsub(
        paste0("wt-\\d"),
        paste0("wt-", i + 1),
        rep_phis_list[[i]]$from, perl = TRUE
      )
      rep_phis_list[[i]]$to <- gsub(
        paste0("wt"),
        paste0("wt-", i),
        rep_phis_list[[i]]$to, perl = TRUE
      )
    }
    rep_phis <- do.call(rbind, rep_phis_list)
    all_phis <- rbind(all_phis, rep_phis)
  }
  # indicate if path has to be bend (for lagged AR effects)
  all_phis$lagged <- ifelse(
    as.numeric(
      gsub("(\\w+)\\((\\d)\\)_(\\d)(\\d)", "\\2", all_phis$Param)
    ) > 1,
    as.numeric(
      gsub("(\\w+)\\((\\d)\\)_(\\d)(\\d)", "\\2", all_phis$Param)
    ) - 1, NA
  )
  # indicate phi type
  all_phis$phi_type <- ifelse(
    grepl("(\\d+)\\1", all_phis$Param),
    "ar", "cl"
  )
  # indicate which paths to draw lagged AR effects
  # (only draw lagged AR effects for first and last time series)
  all_phis$first <- ifelse(
    as.numeric(gsub(".*?(\\d).*", "\\1", all_phis$from)) ==
      min(as.numeric(gsub(".*?(\\d).*", "\\1", all_phis$from))),
    1, 0
  )
  all_phis$last <- ifelse(
    as.numeric(gsub(".*?(\\d).*", "\\1", all_phis$from)) ==
      max(as.numeric(gsub(".*?(\\d).*", "\\1", all_phis$from))),
    1, 0
  )
  # don't draw AR effects of lag > 1 for middle time series
  all_phis$draw <- ifelse(
    all_phis$phi_type == "ar" &
      !is.na(all_phis$lagged) &
      all_phis$first == 0 & all_phis$last == 0,
    0, 1
  )
  # indicate which paths to label
  all_phis$label <- ifelse(
    is.na(all_phis$lagged), 1, ifelse(
      !is.na(all_phis$lagged) & all_phis$phi_type == "ar", 1, 0
    )
  )
  # overwrite names
  all_phis$names <- ifelse(all_phis$label == 1, all_phis$names, "")
  # delete paths that should not be drawn
  all_phis <- all_phis[all_phis$draw == 1, ]
  # loop across all phis
  n_phi <- nrow(all_phis)
  phi_vec <- c()
  for (i in 1:n_phi) {
    phi_vec <- c(
      phi_vec, paste0(
        "\\draw  [path",
        ifelse(
          # decorate with dot on path if parameter is random
          all_phis$isRandom[i] == 1 & all_phis$phi_type[i] == "ar",
          ", postaction = random",
          ifelse(
            all_phis$isRandom[i] == 1 & all_phis$phi_type[i] == "cl",
            ", postaction = random_cl", ""
          )
        ),
        ifelse(
          # set anchor explicitly for lagged AR paths of last time series
          all_phis$phi_type[i] == "ar" &
            !is.na(all_phis$lagged[i]) &
            all_phis$last[i] == 1,
          ", anchor = west", ""
        ),
        ifelse(
          # set anchor explicitly for lagged AR paths of last time series
          all_phis$phi_type[i] == "cl",
          ", sloped, pos = .2", ""
        ), "]  (",
        # select start node conditional on infos$p
        ifelse(infos$p > 1, "eta", "y"), all_phis$from[i],
        # bend path for lagged AR effects
        ifelse(
          !is.na(all_phis$lagged[i]) &
            all_phis$phi_type[i] == "ar" &
            all_phis$first[i] == 1,
          paste0(
            ".north)  |-  +(0em, ",
            all_phis$lagged[i], "em) -|  node  []  {"
          ),
          ifelse(
            !is.na(all_phis$lagged[i]) &
              all_phis$phi_type[i] == "ar" &
              all_phis$last[i] == 1,
            paste0(
              ".south)  |-  +(0em, -",
              all_phis$lagged[i], "em) -|  node  []  {"
            ),
          ")  to node  []  {"
          )
        ), all_phis$names[i], "}  (",
        # select target node conditional on infos$p
        ifelse(
          infos$p > 1, "eta", "y"
        ), all_phis$to[i],
        ifelse(
          !is.na(all_phis$lagged[i]) &
            all_phis$phi_type[i] == "ar" &
            all_phis$first[i] == 1,
          ".north);",
          ifelse(
            !is.na(all_phis$lagged[i]) &
              all_phis$phi_type[i] == "ar" &
              all_phis$last[i] == 1,
            ".south);", ");"
          )
        )
      )
    )
  }
  # paste phis in one string
  phi <- paste(phi_vec, collapse = "\n")

  # do the same for innovation variances
  all_sigmas <- model[grepl("ln.sigma2", model$Param) & grepl("Fix", model$Type), ]
  n_sigma <- nrow(all_sigmas)
  sigma_vec <- c()
  for (i in 1:n_sigma) {
    sigma_vec <- c(
      sigma_vec, paste0(
        "\\draw  [var", ifelse(
          # decorate with dot on path if parameter is random
          all_sigmas[all_sigmas$Param == paste0("ln.sigma2_", i), "isRandom"] == 1,
          ", postaction = random]", "]"
        ),
        "  (delta", i, "t.120)  to node  []  {$\\pi_{", i, "}$}  (delta", i, "t.60);"
      )
    )
  }
  # paste sigmas in one string
  sigma <- paste(sigma_vec, collapse = "\n")
  # do the same for innovation covariances
  all_psis <- model[grepl("ln.sigma_", model$Param) & grepl("Fix", model$Type), ]
  # n_psi <- nrow(all_psis)
  psi_vec <- c()
  if (nrow(all_psis) >= 1) {
    if (infos$q == 1) {
      psi_vec <- paste0(
        "\\draw  [cov", ifelse(
          # decorate with dot on path if parameter is random
          all_psis[all_psis$Param == "ln.sigma_12", "isRandom"] == 1,
          ", postaction = random]", "]"
        ),
        "  (delta1t.0)  to node  []  {$\\pi_{12}$}  (delta2t.0);"
      )
    } else {
      for (i in 1:(infos$q - 1)) {
        for (j in (i + 1):infos$q) {
          psi_vec <- c(
            psi_vec, paste0(
              "\\draw  [cov", ifelse(
                # decorate with dot on path if parameter is random
                all_psis[all_psis$Param == paste0("ln.sigma_", i, j), "isRandom"] == 1,
                ", postaction = random]", "]"
              ),
              "  (delta", i, "t.0)  to node  []  {$\\pi_{", i, j, "}$}  (delta", j, "t.0);"
            )
          )
        }
      }
    }
  }
  # paste psis in one string
  psi <- paste(psi_vec, collapse = "\n")

  # paste together
  within_model <- paste(wm, phi, sigma, psi, sep = "\n")
  # within_model <- paste(wm, sep = "\n")

  # paste together
  within_model <- paste(
    # begin_figure,
    # wm_caption,
    begin_tikz,
    within_model,
    wm_caption,
    end_tikz,
    # end_figure,
    sep = "\n"
  )

  # paste within_model to markdown
  cat(within_model, file = "pathmodel.rmd", append = TRUE)


  # between-model #############################################################

  # between-model caption
  bm_caption <- paste0(
    "% caption\n\\node  [above = 1em, align = center]  at  ",
    "(current bounding box.north)  {Between-model.};"
  )

  # get random effects from within-model fixed effects
  all_bpars <- model[grepl("Fix", model$Type), ]
  # the next line removes non-random effects completely, if desired
  # all_bpars <- model[grepl("Fix", model$Type) & model$isRandom == 1, ]

  # replace dots with nothing
  all_bpars$Param <- gsub("\\.", "", all_bpars$Param)
  all_bpars$names <- gsub(
    # replace underscore digit with latex subscript
    "(\\w+)_(\\d+)", "$\\\\\\1_{\\2}$", gsub(
      # replace lnsigma with log(pi)
      "(lnsigma2|lnsigma)_(\\d+)", "log($\\\\pi_{\\2}$)", gsub(
        # replace uppercase B for between-level latent variables
        "(\\w+)(B)_(\\d)", "$\\\\\\1^{b}_\\3$", gsub(
          # replace underscore digit with latex subscript for phis
          "(\\w+)(\\(\\d\\))_(\\d+)", "$\\\\\\1_{\\3;\\2}$", all_bpars$Param
        )
      )
    )
  )
  # delete rounded brackets from parameters (not allowed in tikz)
  # has to be done after all_bpars$names are generated!
  all_bpars$Param <- gsub("\\(|\\)", "", all_bpars$Param)
  n_bpars <- nrow(all_bpars)
  bpars_vec <- c()
  for (i in 1:n_bpars) {
    if (i == 1) {
      bpars_vec <- c(
        bpars_vec, paste0(
          "\\node  [latent", ifelse(
            # draw in gray if parameter is fixed
            all_bpars[i, "isRandom"] == 0,
            ", color = gray]  (", "]  ("
          ),
          all_bpars[i, "Param"], ")  {", all_bpars$names[i], "};"
        )
      )
    } else {
      bpars_vec <- c(
        bpars_vec, paste0(
          "\\node  [latent", ifelse(
            # draw in gray if parameter is fixed
            all_bpars[i, "isRandom"] == 0,
            ", color = gray]  (", "]  ("
          ),
          all_bpars[i, "Param"], ")  [right = 1em of ",
          all_bpars[i - 1, "Param"],
          "]  {", all_bpars$names[i], "};"
        )
      )
    }
  }
  bpars <- paste(bpars_vec, collapse = "\n")

  # get random effect predictors and predicted outcomes if provided
  if (nrow(infos$RE.PREDS) > 0 | nrow(infos$OUT) > 0) {
    # for predicted random effects
    if (nrow(infos$RE.PREDS) > 0) {
      re_preds <- infos$RE.PREDS[, c(
        "Model", "Level", "Type", "Param", "isRandom", "Lag", "isAR",
        "re_as_dv", "re_preds", "re_no"
      )]
      # rename for easier binding later
      names(re_preds)[8] <- "dv"
      names(re_preds)[9] <- "pred"
      names(re_preds)[10] <- "no"
    } else {
      # if not provided, create empty data frame for binding
      re_preds <- data.frame(matrix(ncol = 10, nrow = 0))
      colnames(re_preds) <- c(
        "Model", "Level", "Type", "Param", "isRandom", "Lag", "isAR",
        "re_as_dv", "re_preds", "re_no"
      )
    }
    # for predicted outcomes
    if (nrow(infos$OUT) > 0) {
      out_preds <- infos$OUT[, c(
        "Model", "Level", "Type", "Param", "isRandom", "Lag", "isAR",
        "Var", "Pred", "out_var_no"
      )]
      # rename for easier binding later
      names(out_preds)[8] <- "dv"
      names(out_preds)[9] <- "pred"
      names(out_preds)[10] <- "no"
    } else {
      # if not provided, create empty data frame for binding
      out_preds <- data.frame(matrix(ncol = 10, nrow = 0))
      colnames(out_preds) <- c(
        "Model", "Level", "Type", "Param", "isRandom", "Lag", "isAR",
        "re_as_dv", "re_preds", "re_no"
      )
    }

    # bind together
    all_preds <- rbind(re_preds, out_preds)
    # replace dots with nothing
    all_preds[, c("Param", "dv", "pred")] <- apply(
      all_preds[, c("Param", "dv", "pred")], 2,
      function(x) gsub("\\.|\\(|\\)", "", x)
    )
    # get additional nodes to draw
    pred_nodes <- c(infos$n_cov_vars, infos$out_var)
    # delete emtpy entries
    pred_nodes <- pred_nodes[!is.na(pred_nodes)]
    # replace special characters
    pred_nodes <- gsub("\\.|\\(|\\)", "", pred_nodes)

    # place coordinates below between-level parameters
    pred_coords_vec <- c()
    for (i in c(1, nrow(all_bpars))) {
      pred_coords_vec <- c(
        pred_coords_vec, paste0(
          "\\node  [coordinate]  (cpred", i, ")  [below = 7.5em of ",
          all_bpars$Param[i], "]  {};"
        )
      )
    }
    pred_coords <- paste(pred_coords_vec, collapse = "\n")

    # place nodes under bpars
    preds_vec <- c()
    for (i in 1:length(pred_nodes)) {
      preds_vec <- c(
        preds_vec, paste0(
          "\\node  [manifest]  (", pred_nodes[i], ")  at  ($(cpred1)!",
          i, "/", length(pred_nodes) + 1, "!(cpred", nrow(all_bpars), ")$)",
          "  {\\textit{", pred_nodes[i], "}};"
        )
      )
    }
    preds <- paste(preds_vec, collapse = "\n")

    # draw paths
    preds_paths_vec <- c()
    n_preds <- nrow(all_preds)
    for (i in 1:n_preds) {
      preds_paths_vec <- c(
        preds_paths_vec, paste0(
          "\\draw  [path]  (", all_preds$pred[i], ")  to node  []  {}  (",
          all_preds$dv[i], ");"
          # ifelse(
          #   # indicate path target anchor
          #   grepl("Outcome", all_preds$Type[i]),
          #   ".north);", ".south);"
          # )
        )
      )
    }
    preds_paths <- paste(preds_paths_vec, collapse = "\n")
  } else {
    pred_coords <- c()
    preds <- c()
    preds_paths <- c()
  }



  bm <- paste(bpars, pred_coords, preds, preds_paths, sep = "\n")

  # draw between-model covariances

  # place coordinates above between-level parameters
  cov_coordinates_vec <- c()
  for (i in 1:n_bpars) {
    cov_coordinates_vec <- c(
      cov_coordinates_vec, paste0(
        "\\node  [coordinate]  (c", i, ")  [above = 2em of ",
        all_bpars[i, "Param"], "]  {};"
      )
    )
  }
  cov_coordinates <- paste(cov_coordinates_vec, collapse = "\n")

  # draw arrows from coordinates to between-level parameters
  cov_arrows_vec <- c()
  for (i in 1:n_bpars) {
    # if parameter is not random, don't draw arrow
    if (all_bpars[i, "isRandom"] == 0) {
      cov_arrows_vec <- c(cov_arrows_vec)
    } else {
      cov_arrows_vec <- c(
        cov_arrows_vec, paste0(
          "\\draw  [path]  (c", i, ")  to node  []  {}  (",
          all_bpars[i, "Param"], ");"
        )
      )
    }
  }
  cov_arrows <- paste(cov_arrows_vec, collapse = "\n")

  # draw line above arrows to connect for covariance
  cov_line <- paste0(
    "\\draw  [thick]  (c",
    # start line above first random between-level parameter
    which(all_bpars[, "isRandom"] == 1)[1],
    ")  --  (c",
    # end line above last random between-level parameter
    which(all_bpars[, "isRandom"] == 1)[length(which(all_bpars[, "isRandom"] == 1))],
    ");"
  )
  # paste together
  covs <- paste(cov_coordinates, cov_arrows, cov_line)
  between_model <- paste(bm, covs, sep = "\n")

  # paste together
  between_model <- paste(
    # begin_figure,
    # bm_caption,
    begin_tikz,
    between_model,
    bm_caption,
    end_tikz,
    # end_figure,
    sep = "\n"
  )

  # paste within_model to markdown
  cat(between_model, file = "pathmodel.rmd", append = TRUE)


  # render markdown input #####################################################

  # render markdown
  rmarkdown::render(
    input = "pathmodel.rmd"
  )

  # optional: store PDF pages as separate pngs ################################
  if (add.png == T) {
    pdftools::pdf_convert(
      pdf = "pathmodel.pdf",
      format = "png",
      pages = NULL,
      filenames = NULL,
      dpi = 720,
      antialias = TRUE,
      opw = "",
      upw = "",
      verbose = FALSE
    )
  }
}
