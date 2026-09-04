# Compare different estimators for a (possibly blocked) cluster randomized trial.

This function calculates the point estimates and SE estimates for a
variety of the estimators implemented by this package.

## Usage

``` r
compare_methods(
  formula,
  data = NULL,
  control_formula = NULL,
  include_MLM = TRUE,
  include_DB = TRUE,
  include_LM = TRUE,
  include_agg = TRUE,
  include_gee = FALSE,
  warn_missing = TRUE,
  include_dumb = FALSE,
  include_disfavored = FALSE,
  patch_data = TRUE,
  handle_singleton_blocks = c("drop", "pool", "fail"),
  include_method_characteristics = TRUE
)
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

- control_formula:

  What variables to control for, in the form of "~ X1 + X2".

- include_MLM:

  Include MLM estimators

- include_DB:

  Include Design-Based estimators (taken from RCTYes documentation and
  prior literature).

- include_LM:

  Include Linear Model-based estimators (including Huber-White SEs,
  etc.)

- include_agg:

  Include estimators applied to aggregated data, aggregated at the
  cluster level.

- include_gee:

  Include GEE estimators. Off by default: the GEE implementation is
  still in development and can fail on rank-deficient designs (see
  \`gee_estimators()\`).

- warn_missing:

  If TRUE, warn when rows are dropped or imputed due to missing data.

- include_dumb:

  Include "dumb" estimators (i.e., those interacted estimators that
  weight by both person and cluster or vice versa).

- include_disfavored:

  Include estimators flagged as "disfavored" in
  \`method_characteristics()\` (odd weighting or known instability).
  Also determines whether unrecognized method names are dropped
  (default) or kept with NA characteristics.

- patch_data:

  If TRUE impute all missing covariates with mean imputation (adding
  dummy variables as needed) via the \`patch_data()\` method. Will drop
  all rows with missing outcome, treatment, or clustering info. If FALSE
  do not do this.

- handle_singleton_blocks:

  How to handle blocks with a treatment or control arm containing only a
  single cluster: "drop" the block, "pool" it with another singleton
  block, or "fail" with an error.

- include_method_characteristics:

  Include details of the methods (target estimands and sampling
  framework assumed) in the return value.

## Value

Dataframe of point estimates and standard errors for each method
considered. If `include_method_characteristics=TRUE` also add some
features of the methods as additional columns.

## Examples

``` r
data( fakeCRT )
compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT )
#>            method  weight   ATE_hat     SE_hat    p_value       df biased
#> 1        AR-FE-db Cluster 0.1852251 0.08740229 0.03537707 189.0000      1
#> 2       AR-FE-het Cluster 0.1852251 0.08735417 0.03527727 189.0000      1
#> 3      AR-FIcw-db Cluster 0.1867940 0.08783710 0.03345371 180.0000      0
#> 4     AR-FIcw-het Cluster 0.1867940 0.08783710 0.03345371 180.0000      0
#> 5      ARpw-FE-db  Person 0.2032865 0.08762329 0.02140859 189.0000      1
#> 6     ARpw-FE-het  Person 0.2032865 0.08660625 0.01994741 189.0000      1
#> 7    ARpw-FIpw-db  Person 0.2044673 0.08567461 0.01700696 180.0000      0
#> 8   ARpw-FIpw-het  Person 0.2044673 0.08600557 0.01743643 180.0000      0
#> 9      LR-FE-crve  Person 0.2032865 0.08660625 0.02015520 157.6603      1
#> 10       LR-FE-db  Person 0.2032865 0.08762329 0.02140859 189.0000      1
#> 11   LR-FIpw-crve  Person 0.2044673 0.08600557 0.01743643 140.2158      0
#> 12     LR-FIpw-db  Person 0.2044673 0.08567461 0.01700696 180.0000      0
#> 13   LRcw-FE-crve Cluster 0.1852251 0.08710924 0.03516777 144.9548      1
#> 14     LRcw-FE-db Cluster 0.1852251 0.08740229 0.03537707 189.0000      1
#> 15 LRcw-FIcw-crve Cluster 0.1867940 0.08780516 0.03338943 134.7822      0
#> 16   LRcw-FIcw-db Cluster 0.1867940 0.08783710 0.03345371 180.0000      0
#> 17         MLM-FE Cluster 0.1931067 0.08674392 0.02720194 186.2165      1
#> 18       MLM-FIcw Cluster 0.1947340 0.08709019 0.02535151 180.0000      1
#> 19         MLM-RE Cluster 0.1816252 0.08666143 0.03744273 187.4170      1
#>    disfavored
#> 1           0
#> 2           0
#> 3           0
#> 4           0
#> 5           0
#> 6           0
#> 7           0
#> 8           0
#> 9           0
#> 10          0
#> 11          0
#> 12          0
#> 13          0
#> 14          0
#> 15          0
#> 16          0
#> 17          0
#> 18          0
#> 19          0
```
