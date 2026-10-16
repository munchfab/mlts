# Create TeX Model Formula from mlts model object

Create TeX Model Formula from mlts model object

## Usage

``` r
mlts_model_formula(
  model,
  file = NULL,
  keep_tex = FALSE,
  ts = NULL,
  covariates = NULL,
  outcomes = NULL
)
```

## Arguments

- model:

  A model built with
  [`mlts_model`](https://munchfab.github.io/mlts/reference/mlts_model.md).

- file:

  An optional string containing the name of the file and file path. Has
  to end with .pdf file format.

- keep_tex:

  Logical. Should the TeX file be kept (additional to the Rmd file)?
  Defaults to `FALSE`.

- ts:

  To be included in future releases. An optional character vector
  containing the names of the time-series variables or indicators.

- covariates:

  To be included in future releases. An optional character vector
  containing the names of the between-level covariates.

- outcomes:

  To be included in future releases. An optional character vector
  containing the names of the between-level outcomes.

## Value

An RMarkdown file that is automatically rendered to a pdf document.

## Examples

``` r
# \donttest{
# build a simple vector-autoregressive mlts model with two time-series variables
var_model <- mlts_model(q = 2)

# create formula from the specified model
mlts_model_formula(model = var_model)
#> 
#> 
#> processing file: formula.rmd
#> 1/3        
#> 2/3 [setup]
#> 3/3        
#> output file: formula.knit.md
#> /opt/hostedtoolcache/pandoc/3.8.3/x64/pandoc +RTS -K512m -RTS formula.knit.md --to latex --from markdown+autolink_bare_uris+tex_math_single_backslash --output formula.tex --lua-filter /home/runner/work/_temp/Library/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /home/runner/work/_temp/Library/rmarkdown/rmarkdown/lua/latex-div.lua --embed-resources --standalone --syntax-highlighting tango --pdf-engine pdflatex --variable graphics --lua-filter /home/runner/work/_temp/Library/rmarkdown/rmarkdown/lua/extract-data-uri.lua --include-in-header /tmp/RtmpIx2bU0/rmarkdown-str1d8b4fe9bd34.html --variable 'geometry:margin=1in' 
#> Warning: error in running command
#> ! sh: 1: pdflatex: not found
#> Error: LaTeX failed to compile formula.tex. See https://yihui.org/tinytex/r/#debugging for debugging tips. See formula.log for more info.
# }
```
