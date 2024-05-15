#' Create Path Diagrams from mlts model object
#'
#' @param model A model built with \code{\link[mlts]{mlts_model}}.
#' @param file An optional string containing the name of the file and file path.
#' Has to end with .pdf file format.
#' @param add_png Logical. Set to `TRUE` to transform created PDF to .png file
#' using `pdftools::pdf_convert`.
#' @param keep_tex Logical. Should the TeX file be kept (additional to the
#' Rmd file)? Defaults to `FALSE`.
#' @param ts To be included in future releases.
#' An optional character vector containing the names of the time-series
#' variables or indicators.
#' @param covariates To be included in future releases.
#' An optional character vector containing the names of the between-level
#' covariates.
#' @param outcomes To be included in future releases.
#' An optional character vector containing the names of the between-level
#' outcomes.
#' @return An RMarkdown file that is automatically rendered to a pdf document.
#' @export
#'
#' @examples
#' \donttest{
#' # build a simple vector-autoregressive mlts model with two time-series variables
#' var_model <- mlts_model(q = 2)
#'
#' # create a pathmodel from the specified model
#' mlts_model_paths(model = var_model)
#' }
mlts_model_paths <- function(model, file = NULL,
                             add_png = FALSE, keep_tex = FALSE,
                             ts = NULL, covariates = NULL,
                             outcomes = NULL) {


  # extract model infos
  infos <- mlts_model_eval(model)

  # extract model data frame
  model <- model

  # get file path, choose top-level directory if not specified
  if (!is.null(file)) {
    rmd_file <- gsub(".pdf", ".rmd", file)
    pdf_file <- file
  } else {
    rmd_file <- "pathmodel.rmd"
    pdf_file <- "pathmodel.pdf"
  }

  # create empty markdown file, delete if already existing
  if (file.exists(rmd_file)) {
    file.remove(rmd_file)
  }
  rmarkdown::draft(file = rmd_file,
                   template = "pathmodel",
                   package = "mlts",
                   edit = FALSE)

  # add keep_tex to YAML frontmatter if desired
  if (keep_tex == TRUE) {
    yaml_keep_tex <- gsub(
      "output: pdf_document", "output: pdf_document:\n\tkeep_tex: true",
      readLines(rmd_file)
    )
    writeLines(yaml_keep_tex, rmd_file)
  }

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
    \\node  [manifest]  (y1t)  {$y_{1,it}$};
    \\node  [latent]  (y1wt)  [above = 2.5em of y1t]  {$y_{1,it}^w$};
    \\node  [latent]  (mu1)  [below = 2.5em of y1t]  {$\\mu_{1,i}$};

    % draw paths
    \\draw  [path]  (y1wt)  to node  []  {}  (y1t);
    \\draw  [path]  (mu1)  to node  []  {}  (y1t);
    "
    if (infos$isLatent == TRUE) { # one latent construct
      # initiate empty vectors to store tikz nodes
      ind_vec <- c() # indicator variables
      # latent variables (within and between)
      wlat <- "\\node  [latent]  (eta1wt)  [above = 4em of c]  {$\\eta^w_{1,it}$};"
      if (all(infos$indicators$btw_factor == 1)) {
        blat <- "\\node  [latent]  (eta1b)  [below = 4em of c]  {$\\eta^b_{1,i}$};"
      } else {
        blat_vec <- c()
        for (i in 1:infos$p) {
          blat_vec[i] <- paste0(
            "\\node  [latent, minimum width = 3em]  (mu1", i,
            ")  [below = 4em of y1", i, "t]  ",
            "{$\\mu_{1", i, ",i}$};"
          )
        }
        blat <- paste(blat_vec, collapse = "\n")
      }
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
            "\\node  [manifest]  (y1", i, "t)  {$y_{1", i, ",it}$};"
          ))
          # first error (labeled)
          epswt_vec <- c(epswt_vec, paste0(
            "\\node  [error]  (eps1", i, "wt)  [above left = 1em of y1",
            i, "t]  {\\scriptsize$\\varepsilon^w_{1", i, ",it}$};"
          ))
          # epsb_vec <- c(epsb_vec, paste0(
          #   "\\node  [error]  (eps", i, "b)  [below left = 1em of y",
          #   i, "t]  {\\scriptsize$\\varepsilon^b_{", i, "}$};"
          # )) # should not be included
        } else {
          # all other indicators
          ind_vec <- c(ind_vec, paste0(
            "\\node  [manifest]  (y1", i, "t)  [right = 1.5em of y1",
            i - 1, "t]  {$y_{1", i, ",it}$};"
          ))
          # residuals
          epswt_vec <- c(epswt_vec, paste0(
            "\\node  [error]  (eps1", i, "wt)  [above left = .75em of y1", i, "t]  {};"
          ))
          if (i == infos$p) {
            # place coordinates between first and last indicator variables
            # of construct
            coord <- paste0(
              "\\node  [coordinate]  (c)  at ($(y11t) !0.5! (y1", i, "t)$)  {};"
            )
            # label last error on between level
            epsb_vec <- c(epsb_vec, ifelse(
              all(infos$indicators$btw_factor == 1),
              paste0(
                "\\node  [error]  (eps1", i, "b)  [below right = 1em of y1",
                i, "t]  {\\scriptsize$\\varepsilon^b_{1", i, ",i}$};"
              ), ""
            ))
          } else {
            epsb_vec <- c(epsb_vec, ifelse(
              all(infos$indicators$btw_factor == 1),
              paste0(
                "\\node  [error]  (eps1", i, "b)  [below right = .75em of y1", i, "t]  {};"
              ), ""
            ))
          }
        }
        # draw paths
        wlat_paths_vec <- c(wlat_paths_vec, paste0(
          "\\draw  [path]  (eta1wt)  to node  [fill = white, anchor = center]
            {\\scriptsize", ifelse(
              # label path with constraint if necessary
              model[model$Param == paste0("lambdaW_1.", i), "Constraint"] == "= 1",
              "$1$}", paste0("$\\lambda^w_{1", i, "}$}")
            ), "  (y1", i, "t.north);"
        ))
        blat_paths_vec <- c(blat_paths_vec, paste0(
          "\\draw  [path]  ", ifelse(
            all(infos$indicators$btw_factor == 1),
            paste0("(eta1b)"), paste0("(mu1", i, ")")
          ), "  to node  ", ifelse(
            all(infos$indicators$btw_factor == 0),
            "[]", "[fill = white, anchor = center]"
          ),  "{", ifelse(
            # label path with constraint if necessary
            model[model$Param == paste0("lambdaB_1.", i), "Constraint"] == "= 1",
            "\\scriptsize $1$", paste0("\\scriptsize $\\lambda^b_{1", i, "}$")
          ), "}  (y1", i, "t.south);"
        ))
        epswt_paths_vec <- c(epswt_paths_vec, paste0(
          "\\draw  [path]  (eps1", i, "wt)  to node  []  {}  (y1", i, "t.north west);"
        ))
        if (i != 1) {
          epsb_paths_vec <- c(epsb_paths_vec, ifelse(
            all(infos$indicators$btw_factor == 1),
            paste0(
              "\\draw  [path]  (eps1", i, "b)  to node  []  {}  (y1", i, "t.south east);"
            ), ""
          ))
        }
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
    \\node  [latent]  (y1wt)  [above = 2.5em of y1t]  {$y_{1,it}^w$};
    \\node  [latent]  (mu1)  [below = 2.5em of y1t]  {$\\mu_{1,i}$};

    % draw nodes
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    infos$q, "}
    \\node  [manifest]  (y\\i t)  [right = 2.5em of y\\lasti t]  {$y_{\\i,it}$};

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    infos$q, "}
    \\node  [latent]  (y\\i wt)  [above = 2.5em of y\\i t]  {$y_{\\i,it}^w$};

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ...,",
    infos$q, "}
    \\node  [latent]  (mu\\i)  [below = 2.5em of y\\i t]  {$\\mu_{\\i,i}$};

    % draw paths
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
    infos$q, "}
    \\draw  [path]  (y\\i wt)  to node  []  {}  (y\\i t);

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
    infos$q, "}
    \\draw  [path]  (mu\\i)  to node  []  {}  (y\\i t);"
    )
    if (infos$isLatent == TRUE) {
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
          "]  {$\\eta^w_{", i, ",it}$};"
        ))
        if (all(infos$indicators[infos$indicators$q == i, "btw_factor"] == 1)) {
          # if common between-factor is modeled, use eta
          blat_vec <- c(blat_vec, paste0(
            "\\node  [latent]  (etab", i, ")  [below = 4em of c", i,
            "]  {$\\eta^b_{", i, ",i}$};"
          ))
        } else {
          # if no between-factor is modeled, use mu
          for (j in 1:infos$p[i]) {
            blat_vec <- c(blat_vec, paste0(
              "\\node  [latent, minimum width = 3em]  (mu", i, j,
              ")  [below = 3em of y", i, j,
              "t]  {$\\mu_{", i, j, ",i}$};"
            ))
          }
        }
        # if one construct only has 1 indicator, place coordinate at indicator node
        if (infos$p[i] == 1) {
          coord_vec <- c(coord_vec, paste0(
            "\\node  [coordinate]  (c", i,
            ")  at (y", i, "1t)  {};"
          ))
        }
        for (j in 1:infos$p[i]) {
          if (i == 1 & j == 1) {
            # first indicator variable
            ind_vec <- c(ind_vec, paste0(
              "\\node  [manifest]  (y", i, j, "t)  {$y_{", i, j, ",it}$};"
            ))
            # first error (labeled)
            if (infos$p[i] == 1) {
              epswt_vec <- epswt_vec
              epsb_vec <- epsb_vec
              epswt_paths_vec <- epswt_paths_vec
              epsb_paths_vec <- epsb_paths_vec
            } else {
              epswt_vec <- c(epswt_vec, paste0(
                "\\node  [error]  (eps", i, j, "wt)  [above left = 1em of y",
                i, j, "t]  {\\scriptsize$\\varepsilon^w_{", i, j, ",it}$};"
              ))
              # epsb_vec <- c(epsb_vec, paste0(
              #   "\\node  [error]  (eps", i, j, "b)  [below right = 1em of y",
              #   i, j, "t]  {\\scriptsize$\\varepsilon^b_{", i, j, "}$};"
              # )) # should not be included
              epswt_paths_vec <- c(epswt_paths_vec, paste0(
                "\\draw  [path]  (eps", i, j, "wt)  to node  []  {}  (y", i, j, "t.north west);"
              ))
              # epsb_paths_vec <- c(epsb_paths_vec, paste0(
              #   "\\draw  [path]  (eps", i, j, "b)  to node  []  {}  (y", i, j, "t.south east);"
              # ))
            }
          } else if (i > 1 & j == 1) {
            # first indicator variable of next construct
            ind_vec <- c(ind_vec, paste0(
              "\\node  [manifest]  (y", i, j, "t)  [right = 2.5em of y",
              i - 1, infos$p[i - 1], "t]  {$y_{", i, j, ",it}$};"
            ))
            if (infos$p[i] == 1) {
              epswt_vec <- epswt_vec
              epsb_vec <- epsb_vec
              epswt_paths_vec <- epswt_paths_vec
              epsb_paths_vec <- epsb_paths_vec
            } else {
              epswt_vec <- c(epswt_vec, paste0(
                "\\node  [error]  (eps", i, j, "wt)  [above left = .75em of y", i, j, "t]  {};"
              ))
              # epsb_vec <- c(epsb_vec, paste0(
              #   "\\node  [error]  (eps", i, j, "b)  [below left = .75em of y", i, j, "t]  {};"
              # ))
              epswt_paths_vec <- c(epswt_paths_vec, paste0(
                "\\draw  [path]  (eps", i, j, "wt)  to node  []  {}  (y", i, j, "t.north west);"
              ))
              # epsb_paths_vec <- c(epsb_paths_vec, paste0(
              #   "\\draw  [path]  (eps", i, j, "b)  to node  []  {}  (y", i, j, "t.south west);"
              # ))
            }
          } else if (i == infos$q & j == infos$p[i]) {
            epswt_vec <- c(epswt_vec, paste0(
              "\\node  [error]  (eps", i, j, "wt)  [above left = .75em of y", i, j, "t]  {};"
            ))
            epswt_paths_vec <- c(epswt_paths_vec, paste0(
              "\\draw  [path]  (eps", i, j, "wt)  to node  []  {}  (y", i, j, "t.north west);"
            ))
            epsb_vec <- c(epsb_vec, ifelse(
              all(infos$indicators[infos$indicators$q == i, "btw_factor"] == 1),
              paste0(
                "\\node  [error]  (eps", i, j, "b)  [below right = 1em of y",
                i, j, "t]  {\\scriptsize$\\varepsilon^b_{", i, j, ",i}$};"
              ), ""
            ))
            epsb_paths_vec <- c(epsb_paths_vec, ifelse(
              all(infos$indicators[infos$indicators$q == i, "btw_factor"] == 1),
              paste0(
                "\\draw  [path]  (eps", i, j, "b)  to node  []  {}  (y", i, j, "t.south east);"
              ), ""
            ))
            # last indicator of last construct
            ind_vec <- c(ind_vec, paste0(
              "\\node  [manifest]  (y", i, j, "t)  [right = 1.5em of y",
              i, j - 1, "t]  {$y_{", i, j, ",it}$};"
            ))
            # place coordinates between first and last indicator variables
            # of construct i
            coord_vec <- c(coord_vec, paste0(
              "\\node  [coordinate]  (c", i,
              ")  at ($(y", i, "1t) !0.5! (y", i, j, "t)$)  {};"
            ))
          } else {
            # all other indicators
            ind_vec <- c(ind_vec, paste0(
              "\\node  [manifest]  (y", i, j, "t)  [right = 1.5em of y",
              i, j - 1, "t]  {$y_{", i, j, ",it}$};"
            ))
            epswt_vec <- c(epswt_vec, paste0(
              "\\node  [error]  (eps", i, j, "wt)  [above left = .75em of y", i, j, "t]  {};"
            ))
            epsb_vec <- c(epsb_vec, ifelse(
              all(infos$indicators[infos$indicators$q == i, "btw_factor"] == 1),
              paste0(
                "\\node  [error]  (eps", i, j, "b)  [below right = .75em of y", i, j, "t]  {};"
              ), ""
            ))
            if (j == infos$p[i]) {
              # place coordinates between first and last indicator variables
              # of construct i
              coord_vec <- c(coord_vec, paste0(
                "\\node  [coordinate]  (c", i,
                ")  at ($(y", i, "1t) !0.5! (y", i, j, "t)$)  {};"
              ))
            }
            epswt_paths_vec <- c(epswt_paths_vec, paste0(
              "\\draw  [path]  (eps", i, j, "wt)  to node  []  {}  (y", i, j, "t.north west);"
            ))
            epsb_paths_vec <- c(epsb_paths_vec, ifelse(
              all(infos$indicators[infos$indicators$q == i, "btw_factor"] == 1),
              paste0(
                "\\draw  [path]  (eps", i, j, "b)  to node  []  {}  (y",
                i, j, "t.south east);"
              ), ""
            ))
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
          if (all(infos$indicators[infos$indicators$q == i, "btw_factor"] == 1)) {
            blat_paths_vec <- c(blat_paths_vec, paste0(
              "\\draw  [path]  (etab",
              i, ")  to node  [fill = white, anchor = center]  {\\scriptsize",
              ifelse(
                # label path with constraint if necessary
                model[model$Param == paste0("lambdaB_", i, ".", j), "Constraint"] == "= 1",
                "$1$}", paste0("$\\lambda^b_{", i, j, "}$}")
              ), "  (y", i, j, "t.south);"
            ))
          } else {
            blat_paths_vec <- c(blat_paths_vec, paste0(
              "\\draw  [path]  (mu",
              i, j, ")  to node  []  {}  (y", i, j, "t.south);"
            ))
          }
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
  cat(decomposition, file = rmd_file, append = TRUE)


  # within-model ##############################################################

  # within-model caption
  wm_caption <- paste0(
    "% caption\n\\node  [above = 1em, align = center]  at  ",
    "(current bounding box.north)  {Within-level model.};"
  )

  # draw nodes for q time-series constructs
  if (infos$q == 1) {
    wm <- paste0(
      "% draw within-level structural model
      \\node  [latent]  (y1wt-1)  {$y_{1,i(t-1)}^w$};
      \\node  [latent]  (y1wt)  [right = 5em of y1wt-1]  {$y_{1,it}^w$};
      \\node  [latent]  (zeta1t)  [right = 2.5em of y1wt]  {$\\zeta_{1,it}$};",

      ifelse(
        infos$maxLag != 1,
        paste0(
          "\\foreach \\j [remember = \\j as \\lastj (initially 1)] in {2, ..., ",
          infos$maxLag, "}
        \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ..., ",
          infos$q, "}
        \\node  [latent]  (y\\i wt-\\j)  [left = 5em of y\\i wt-\\lastj]  {$y_{\\i ,i(t-\\j)}^w$};"
        ),
        paste0("")
      ),

      "% draw paths
      \\draw  [path]  (zeta1t)  to node  []  {}  (y1wt);
    "
    )
    if (infos$isLatent == TRUE) {
      wm <- paste0(
        "% draw within-level structural model
        \\node  [latent]  (eta1wt-1)  {$\\eta_{1,i(t-1)}^w$};
        \\node  [latent]  (eta1wt)  [right = 5em of eta1wt-1]  {$\\eta_{1,it}^w$};
        \\node  [latent]  (zeta1t)  [right = 2.5em of eta1wt]  {$\\zeta_{1,it}$};",

        ifelse(
        infos$maxLag != 1,
          paste0(
            "\\foreach \\j [remember = \\j as \\lastj (initially 1)] in {2, ..., ",
            infos$maxLag, "}
          \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ..., ",
            infos$q, "}
          \\node  [latent]  (eta\\i wt-\\j)  [left = 5em of eta\\i wt-\\lastj]  {$\\eta_{\\i ,i(t-\\j)}^w$};"
          ),
          paste0("")
        ),

        "% draw paths
        \\draw  [path]  (zeta1t)  to node  []  {}  (eta1wt);
      "
      )
    }
  } else { # for q > 1
    wm <- paste0(
    "
    % draw within-level structural model
    \\node  [latent]  (y1wt-1)  {$y_{1,i(t-1)}^w$};
    \\node  [latent]  (y1wt)  [right =5em of y1wt-1]  {$y_{1,it}^w$};
    \\node  [latent]  (zeta1t)  [right =2.5em of y1wt]  {$\\zeta_{1,it}$};

    % draw nodes
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
      infos$q, "}
    \\node  [latent]  (y\\i wt-1)  [below = 5em of y\\lasti wt-1]  {$y_{\\i ,i(t-1)}^w$};",

    ifelse(
      infos$maxLag != 1,
      paste0(
        "\\foreach \\j [remember = \\j as \\lastj (initially 1)] in {2, ..., ",
          infos$maxLag, "}
        \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ..., ",
          infos$q, "}
        \\node  [latent]  (y\\i wt-\\j)  [left = 5em of y\\i wt-\\lastj]  {$y_{\\i ,i(t-\\j)}^w$};"
      ),
      paste0("")
    ),

    "\\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    infos$q, "}
    \\node  [latent]  (y\\i wt)  [below = 5em of y\\lasti wt]  {$y_{\\i ,it}^w$};

    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
    infos$q, "}
    \\node  [latent]  (zeta\\i t)  [below = 5em of zeta\\lasti t]  {$\\zeta_{\\i ,it}$};

    % draw paths from residuals to yiwt
    \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
    infos$q, "}
    \\draw  [path]  (zeta\\i t)  to node  []  {}  (y\\i wt);
    ")
    if (infos$isLatent == TRUE) {
      wm <- paste0(
      "
      % draw within-level structural model
      \\node  [latent]  (eta1wt-1)  {$\\eta_{1,i(t-1)}^w$};
      \\node  [latent]  (eta1wt)  [right = 5em of eta1wt-1]  {$\\eta_{1,it}^w$};
      \\node  [latent]  (zeta1t)  [right = 2.5em of eta1wt]  {$\\zeta_{1,it}$};

      % draw nodes
      \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
        infos$q, "}
      \\node  [latent]  (eta\\i wt-1)  [below = 5em of eta\\lasti wt-1]  {$\\eta_{\\i ,i(t-1)}^w$};",

      ifelse(
        infos$maxLag != 1,
        paste0(
          "\\foreach \\j [remember = \\j as \\lastj (initially 1)] in {2, ..., ",
            infos$maxLag, "}
          \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ..., ",
            infos$q, "}
          \\node  [latent]  (eta\\i wt-\\j)  [left = 5em of eta\\i wt-\\lastj]  {$\\eta_{\\i ,i(t-\\j)}^w$};"
        ),
        paste0("")
      ),

      "\\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
        infos$q, "}
      \\node  [latent]  (eta\\i wt)  [below = 5em of eta\\lasti wt]  {$\\eta_{\\i ,it}^w$};

      \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {2, ..., ",
        infos$q, "}
      \\node  [latent]  (zeta\\i t)  [below = 5em of zeta\\lasti t]  {$\\zeta_{\\i ,it}$};

      % draw paths from residuals to etaiwt
      \\foreach \\i [remember = \\i as \\lasti (initially 1)] in {1, ...,",
        infos$q, "}
      \\draw  [path]  (zeta\\i t)  to node  []  {}  (eta\\i wt);
      ")
    }
  }

  # store number of phi-parameters for loops
  all_phis <- model[grepl("phi", model$Param) & grepl("Fix", model$Type), ]
  all_phis$names <- ifelse(
    all_phis$isRandom == 1,
    # add i index if parameter is random
    gsub(
      # replace underscore digit with latex subscript
      "(\\w+)(\\(\\d\\))_(\\d+)", "$\\\\\\1_{\\2\\3,i}$", all_phis$Param
    ),
    # else don't
    gsub(
      # replace underscore digit with latex subscript
      "(\\w+)(\\(\\d\\))_(\\d+)", "$\\\\\\1_{\\2\\3}$", all_phis$Param
    )
  )
  # generate path start and end
  all_phis$from <- gsub(
    "(\\w+)\\((\\d)\\)_(\\d)(\\d)", "\\4wt-\\2", all_phis$Param
  )
  all_phis$to <- gsub(
    "(\\w+)\\((\\d)\\)_(\\d)(\\d)", "\\3wt", all_phis$Param
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
        # select start node conditional on isLatent
        ifelse(infos$isLatent == TRUE, "eta", "y"), all_phis$from[i],
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
        # select target node conditional on isLatent
        ifelse(
          infos$isLatent == TRUE, "eta", "y"
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
  # all_sigmas <- model[grepl("ln.sigma2|sigma_", model$Param) & grepl("Fix", model$Type), ]
  all_sigmas <- model[
    startsWith(model$Param, "ln.sigma2") | startsWith(model$Param, "sigma_") &
      grepl("Fix", model$Type),
  ]
  n_sigma <- nrow(all_sigmas)
  sigma_vec <- c()
  for (i in 1:n_sigma) {
    sigma_vec <- c(
      sigma_vec, paste0(
        "\\draw  [var", ifelse(
          # decorate with dot on path if parameter is random
          # all_sigmas[all_sigmas$Param == paste0("ln.sigma2_", i), "isRandom"] == 1,
          all_sigmas[i, "isRandom"] == 1,
          ", postaction = random]", "]"
        ),
        "  (zeta", i, "t.30)  to node  []  {$\\sigma", ifelse(
          # if ln.sigma present, use variance, else standard deviation (?)
          grepl("ln", all_sigmas$Param[i]),
          "^2", ""
        ),
        "_{\\zeta_{", i, "}", ifelse(
          # add i index if parameter is random
          all_sigmas$isRandom[i] == 1, ",i", ""
        ), "}$}  (zeta", i, "t.330);"
      )
    )
  }
  # paste sigmas in one string
  sigma <- paste(sigma_vec, collapse = "\n")
  # do the same for innovation covariances
  # all_psis <- model[grepl("ln.sigma_|r.zeta", model$Param) & grepl("Fix", model$Type), ]
  all_psis <- model[
    startsWith(model$Param, "ln.sigma_") | startsWith(model$Param, "r.zeta") &
      grepl("Fix", model$Type),
  ]
  # n_psi <- nrow(all_psis)
  psi_vec <- c()
  # if (nrow(all_psis) >= 1) {
  #   if (infos$q == 1) { # is this needed?
  #     psi_vec <- paste0(
  #       "\\draw  [cov", ifelse(
  #         # decorate with dot on path if parameter is random
  #         all_psis[all_psis$Param == "ln.sigma_12", "isRandom"] == 1,
  #         ", postaction = random]", "]"
  #       ),
  #       "  (zeta1t.0)  to node  []  {$\\psi_{12}$}  (zeta2t.0);"
  #     )
  #   } else {
  #     for (i in 1:(infos$q - 1)) {
  #       for (j in (i + 1):infos$q) {
  #         if (any(grepl("ln.sigma", all_psis$Param)) == TRUE) {
  #           psi_vec <- c(
  #             psi_vec, paste0(
  #               "\\draw  [cov", ifelse(
  #                 # decorate with dot on path if parameter is random
  #                 all_psis[all_psis$Param == paste0("ln.sigma_", i, j), "isRandom"] == 1,
  #                 ", postaction = random]", "]"
  #               ),
  #               "  (zeta", i, "t.0)  to node  []  {$\\psi_{", i, j, "}$}  (zeta", j, "t.0);"
  #             )
  #           )
  #         } else {
  #           psi_vec <- c(
  #             psi_vec, paste0(
  #               "\\draw  [cov]",
  #               "  (zeta", i, "t.0)  to node  []  {$\\psi_{", i, j,
  #               "}$}  (zeta", j, "t.0);"
  #             )
  #           )
  #         }
  #       }
  #     }
  #   }
  # }
  if (nrow(all_psis) >= 1) {
    for (i in 1:(infos$q - 1)) {
      for (j in (i + 1):infos$q) {
        # if random innovation covariance is specified (only for q == 2),
        # draw additional factor for innovation covariance
        if (any(grepl("ln.sigma", all_psis$Param)) == TRUE) {
          psi_vec <- c(
            psi_vec,
            paste0(
              "\\node  [coordinate]  (cetainno", i, j,
              ")  at  ($(zeta", i, "t) !0.5! (zeta", j, "t)$)  {};"
            ),
            paste0(
              "\\node  [latent]  (etainno", i, j,
              ")  at  (cetainno", i, j,
              ")  {$\\eta_{\\zeta_{", i, j, "}, it}$};"
            ),
            paste0(
              "\\draw  [path]",
              "  (etainno", i, j, ")  to node  ",
              "[fill = white, anchor = center]  {\\scriptsize$",
              infos$inno_cov_load[i],
              "$}  (", ifelse(infos$isLatent == TRUE, "eta", "y"),
              i, "wt);"
            ),
            paste0(
              "\\draw  [path]",
              "  (etainno", i, j, ")  to node  ",
              "[fill = white, anchor = center]  {\\scriptsize$",
              infos$inno_cov_load[j],
              "$}  (", ifelse(infos$isLatent == TRUE, "eta", "y"),
              j, "wt);"
            ),
            paste0(
              "\\draw  [var, postaction = random]",
              "  (etainno", i, j, ".30)  to node  []  {$\\sigma",
              "_{\\zeta_{", i, j, "},i}$}  (etainno", i, j, ".330);"
            )
          )
        } else {
          psi_vec <- c(
            psi_vec, paste0(
              "\\draw  [cov]",
              "  (zeta", i, "t.320)  to node  []  {$\\sigma_{\\zeta_{", i, j,
              "}}$}  (zeta", j, "t.40);"
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
  cat(within_model, file = rmd_file, append = TRUE)


  # between-model #############################################################

  # between-model caption
  bm_caption <- paste0(
    "% caption\n\\node  [above = 1em, align = center]  at  ",
    "(current bounding box.north)  {Between-level model.};"
  )

  # get random effects from within-model fixed effects
  all_bpars <- model[grepl("Fix", model$Type), ]
  # the next line removes non-random effects completely, if desired
  # all_bpars <- model[grepl("Fix", model$Type) & model$isRandom == 1, ]
  # replace dots with nothing
  all_bpars$Param <- gsub("\\.", "", all_bpars$Param)
  # replace rzeta with psi (ugly fix but ok)
  all_bpars$Param <- gsub("rzeta", "sigma", all_bpars$Param)
  all_bpars$names <- gsub(
    # replace underscore digit with latex subscript
    "(\\w+)_(\\d+)", "$\\\\\\1_{\\2}$", gsub(
      # replace lnsigma with ln(sigma^2)
      "(lnsigma2)_(\\d+)", "ln($\\\\sigma^2_{\\\\zeta_{\\2}}$)", gsub(
        # replace uppercase B for between-level latent variables
        "(\\w+)(B)_(\\d)", "$\\\\\\1^{b}_{\\3}$", gsub(
          # replace underscore digit with latex subscript for phis
          "(\\w+)(\\(\\d\\))_(\\d+)", "$\\\\\\1_{\\2\\3}$", gsub(
            # replace rzeta with psi in case of fixed
            # innovation covariance (ugly fix but ok)
            "(sigma)_(\\d+)", "$\\\\sigma_{\\\\zeta_{\\2}}$", gsub(
              # replace lnsigma12 (random innovation covariance)
              # with ln(psi)
              "(lnsigma)_(\\d+)", "ln($\\\\sigma_{\\\\zeta_{\\2}}$)", all_bpars$Param
            )
          )
        )
      )
    )
  )
  # add i index to names if parameter is random
  all_bpars$names <- ifelse(
    all_bpars$isRandom == 1,
    gsub("}\\$", ",i}\\$", all_bpars$names),
    all_bpars$names
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
    pred_nodes <- c(infos$n_cov_vars, infos$out_var, infos$n_z_vars)
    # delete duplicate nodes
    pred_nodes <- pred_nodes[!duplicated(pred_nodes) == TRUE]
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
  cat(between_model, file = rmd_file, append = TRUE)


  # render markdown input #####################################################

  # render markdown
  rmarkdown::render(
    input = rmd_file
  )

  # optional: store PDF pages as separate pngs ################################
  if (add_png == T) {
    pdftools::pdf_convert(
      pdf = pdf_file,
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
