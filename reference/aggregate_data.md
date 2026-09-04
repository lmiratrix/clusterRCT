# Aggregate data to the cluster level.

Also implements checks to ensure that all clusters are entirely treated
or not treated.

## Usage

``` r
aggregate_data(data, formula = NULL, control_formula = NULL)
```

## Arguments

- data:

  Dataframe with 'clusterID' and 'Z' as columns. 'blockID' optional
  column. This is data in the "canonical form".

- formula:

  Notation for Y ~ Z \| clusterID \| blockID. If NULL, \`data\` is
  assumed to already be in canonical form.

- control_formula:

  What variables to control for, in the form of "~ X1 + X2". These will
  be averaged (or dummy-expanded and averaged, for categorical
  covariates) to the cluster level.

## Value

tibble of cluster-aggregated data, including Ybar, n, blockID,
clusterID, and Z

## Details

If covariates specified, categorical covariates will be converted to
dummy variables and then averaged (even if they are only level-2 or
level-3 covariates)

Will convert data to canonical form with Yobs, Z, clusterID and blockID
as the names of outcome, treatment, cluster id and block id.
