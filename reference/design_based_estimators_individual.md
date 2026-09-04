# Estimate ATE via design-based approach

This function follows \`design_based_estimators\`, but runs the
regression on the individual level data.

## Usage

``` r
design_based_estimators_individual(
  formula,
  data = NULL,
  control_formula = NULL,
  weight = c("Person", "Cluster"),
  include_block_estimates = FALSE
)
```

## Arguments

- formula:

  Formula for outcome and treatment and nesting. If NULL, data is
  assumed to be in canonical form (see vignette for further discussion).

- data:

  Data frame holding the outcome, treatment, and clustering/blocking
  columns referenced by \`formula\` (or already in canonical form if
  \`formula\` is NULL).

- control_formula:

  What variables to control for, in the form of "~ X1 + X2".

- weight:

  "Person" to weight the regression by individual (unweighted OLS),
  "Cluster" to weight so each cluster contributes equally (inverse
  cluster-size weights).

- include_block_estimates:

  If TRUE, also return the block-by-block estimates used to build the
  overall estimate, in addition to the aggregate estimate.

## Value

tibble of estimates using different varieties of the methods described
in the paper.

## Details

Degrees of freedom: unlike
[`design_based_estimators()`](https://lmiratrix.github.io/clusterRCT/reference/design_based_estimators.md),
this function works directly on individual-level data, so
`control_formula` may mix level-1 (individual) and level-2 (cluster)
covariates. The degrees of freedom adjustment here uses
`number_level2_controls()`, which counts only the level-2 covariates –
level-1 covariates do not cost the same degrees of freedom in this
design-based framework (see the technical supplement's benchmark degrees
of freedom of `J - K - g - 1`, where `g` is the number of level-2
covariates beyond the treatment indicator).
