# Given data, return a data frame with the structure of the data

Given data, return a data frame with the structure of the data

## Usage

``` r
get_structure(formula = NULL, data)
```

## Arguments

- formula:

  Formula of the form \`Yobs ~ Z \| clusterID\` (or \`Yobs ~ Z \|
  clusterID \| blockID\` if blocked). If NULL, \`data\` is assumed to
  already be in canonical form (columns \`Yobs\`, \`Z\`, \`clusterID\`,
  and optionally \`blockID\`).

- data:

  Data frame holding the outcome, treatment, and clustering/blocking
  columns referenced by \`formula\` (or already in canonical form if
  \`formula\` is NULL).

## Value

Dataframe listing blockID, clusterID, treatment status, and number of
units in the cluster.
