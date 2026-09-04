# Calculate number of level 1 and level 2 variables

Given a control formula, determine which variables are level 1 and which
are level 2. Level 1 variables vary within clusters, level 2 do not. If
formula and data are provided, will first convert to canonical form. If
just_count is TRUE, will return a named vector with counts of level 1
and level 2 variables.

## Usage

``` r
count_covariates(
  data,
  control_formula,
  cluster_formula = NULL,
  just_count = FALSE,
  pure = TRUE
)
```

## Arguments

- data:

  Data frame with outcome, treatment, clusterID, blockID, and
  covariates.

- control_formula:

  Formula with covariates only (no outcome or treatment).

- cluster_formula:

  Formula for cluster ID, of form ~ CID or ~ CID \| BID. In the latter
  case, will make interaction to ensure all clusters have unique ID. If
  null, assumes a clusterID and blockID column exist in the data.

- just_count:

  If TRUE, return just counts of level 1 and level 2 variables.

- pure:

  TRUE means count level 2 covariates ONLY if they are ONLY level 2.
  FALSE means covariates that are both level 1 and level 2 get counted
  twice, once for level 1 and once for level 2.

## Value

If just_count is FALSE (default), a tibble with columns:

- variable:

  Variable name.

- level:

  Either "Level-1" or "Level-2".

If just_count is TRUE, a named vector with counts of level 1 and level 2
variables.
