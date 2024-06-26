% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mlts_model_formula.R
\name{mlts_model_formula}
\alias{mlts_model_formula}
\title{Create TeX Model Formula from mlts model object}
\usage{
mlts_model_formula(
  model,
  file = NULL,
  keep_tex = FALSE,
  ts = NULL,
  covariates = NULL,
  outcomes = NULL
)
}
\arguments{
\item{model}{A model built with \code{\link[mlts]{mlts_model}}.}

\item{file}{An optional string containing the name of the file and file path.
Has to end with .pdf file format.}

\item{keep_tex}{Logical. Should the TeX file be kept (additional to the
Rmd file)? Defaults to \code{FALSE}.}

\item{ts}{To be included in future releases.
An optional character vector containing the names of the time-series
variables or indicators.}

\item{covariates}{To be included in future releases.
An optional character vector containing the names of the between-level
covariates.}

\item{outcomes}{To be included in future releases.
An optional character vector containing the names of the between-level
outcomes.}
}
\value{
An RMarkdown file that is automatically rendered to a pdf document.
}
\description{
Create TeX Model Formula from mlts model object
}
\examples{
\donttest{
# build a simple vector-autoregressive mlts model with two time-series variables
var_model <- mlts_model(q = 2)

# create formula from the specified model
mlts_model_formula(model = var_model)
}

}
