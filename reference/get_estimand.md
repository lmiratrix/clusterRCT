# Get the estimand for a given method.

Get the estimand for a given method.

## Usage

``` r
get_estimand(method, simple = TRUE)
```

## Arguments

- method:

  Character vector of method name(s), as they appear in
  \`method_characteristics()\`'s \`method\` column.

- simple:

  If TRUE, collapse "Cluster-Block"/"Person-Block" weight labels down to
  "Cluster"/"Person". If FALSE, collapse only
  "Cluster-Cluster"/"Person-Person" (leaving "-Block" labels distinct).

## See also

method_characteristics()
