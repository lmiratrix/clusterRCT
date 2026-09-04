# Get table of characteristics of all the methods implemented in this package.

Get table of characteristics of all the methods implemented in this
package.

## Usage

``` r
method_characteristics(include_weight = TRUE)
```

## Arguments

- include_weight:

  If TRUE (default), include the \`weight\` column (the estimand each
  method targets). If FALSE, drop it.

## Value

Dataframe with columns for method, weight, population, biased, and
disfavored (whether we do not like the estimator due to odd weighting or
instability).

## Examples

``` r
method_characteristics()
#> # A tibble: 52 × 5
#>    method    weight  biased blocked disfavored
#>    <chr>     <chr>    <dbl>   <dbl>      <dbl>
#>  1 AR-db     Cluster      0       0          0
#>  2 AR-het    Cluster      0       0          0
#>  3 ARpw-db   Person       0       0          0
#>  4 ARpw-het  Person       0       0          0
#>  5 LR-crve   Person       0       0          0
#>  6 LR-db     Person       0       0          0
#>  7 LRcw-crve Cluster      0       0          0
#>  8 LRcw-db   Cluster      0       0          0
#>  9 MLM       Cluster      1       0          0
#> 10 DB_HT     Cluster      0       1          1
#> # ℹ 42 more rows
```
