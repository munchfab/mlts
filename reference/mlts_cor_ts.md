# Correlations for time series variables

Computes correlations between time series variables at the
between-person or within-person level.

## Usage

``` r
mlts_cor_ts(
  data,
  ts,
  id,
  type = "between",
  as_mat = TRUE,
  within_agg_fn = NULL
)
```

## Arguments

- data:

  A data.frame containing the time series data.

- ts:

  Character vector of column names representing the time series
  variables. Must contain at least two variables.

- id:

  Character string specifying the subject identifier column.

- type:

  Character string specifying the level of analysis. Must be either
  `"between"` or `"within"`.

- as_mat:

  Logical. If `TRUE`, returns correlation matrices. If `FALSE`, returns
  a data.frame with pairwise correlations and confidence intervals.

- within_agg_fn:

  Optional function used to aggregate within-person correlations across
  individuals (e.g., `mean`, `median`). Only used when `type = "within"`
  and `as_mat = TRUE`.

## Value

If `type = "between"`:

- If `as_mat = TRUE`: A correlation matrix.

- If `as_mat = FALSE`: A data.frame with variables, correlations, and
  confidence intervals.

If `type = "within"`:

- If `as_mat = TRUE` and `within_agg_fn = NULL`: A list of individual
  correlation matrices.

- If `as_mat = TRUE` and `within_agg_fn` is provided: A single
  aggregated correlation matrix.

- If `as_mat = FALSE`: A data.frame with individual-level correlations
  and confidence intervals.

## Details

At the between-person level, correlations are computed based on person
means. At the within-person level, correlations are computed within
clusters and can optionally be aggregated across clusters using a
summary function.

Correlations are computed using Pearson's correlation coefficient with
pairwise complete observations. Confidence intervals are obtained using
[`cor.test`](https://rdrr.io/r/stats/cor.test.html). If estimation
fails, `NA` values are returned.

## Examples

``` r
if (FALSE) { # \dontrun{
# Specify a three-indicator AR(1) model
model <- mlts_model(q = 1, p = 3)

# Generate data under the specified model
df <- mlts_sim(model = model, N = 50, TP = 100, default = T)$data

# Between-person correlations
mlts_cor_ts(data = df, ts = c("Y1.1", "Y1.2", "Y1.3"), id = "ID")

# Within-person correlations (list of matrices)
mlts_cor_ts(df, ts = c("Y1.1", "Y1.2", "Y1.3"), id = "ID", type = "within")

# Aggregated within-person correlations
mlts_cor_ts(df, ts = c("Y1.1", "Y1.2", "Y1.3"), id = "ID",
            type = "within", within_agg_fn = mean)
} # }
```
