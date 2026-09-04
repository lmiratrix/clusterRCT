# Summarize set of compare_methods results across estimators

Summarize set of compare_methods results across estimators

## Usage

``` r
summarize_simulation_results(rps, summarize_results = "method")
```

## Arguments

- rps:

  A list of data frames, each containing the results of a
  compare_methods call.

- summarize_results:

  If "cross", summarize the results by calculating range statistics
  across the set, grouping estimators by estimand within each runID. If
  "cross-agg", aggregate the sets of ranges after doing this and return
  summary of that. If "method", summarize the results for each method.
  Otherwise just stack the initial set of results (if not already
  stacked)

## Value

A data frame containing the summarized results.
