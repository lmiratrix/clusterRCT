# Estimate ATE via design-based approach

This follows "Design-Based Ratio Estimators and Central Limit Theorems
for Clustered, Blocked RCTs" by Peter Z. Schochet, Nicole E. Pashley,
Luke W. Miratrix, and Tim Kautz.

## Usage

``` r
design_based_estimators(
  formula,
  data = NULL,
  control_formula = NULL,
  weight = c("Person", "Cluster"),
  aggregated = FALSE,
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

- aggregated:

  TRUE means data is already aggregated (and in canonical form). FALSE
  means it is not.

- include_block_estimates:

  If TRUE, also return the block-by-block estimates used to build the
  overall estimate, in addition to the aggregate estimate.

## Value

tibble of estimates using different varieties of the methods described
in the paper.

## Details

It runs a weighted linear regression, as discussed in that paper. The
paper discusses two possibilities: running on the aggregate, or running
on the individual.

In this implementation, we run the regression on the cluster-aggregated
data.

Degrees of freedom: this function receives already cluster-aggregated
data, so every covariate in `control_formula` is, by construction, a
cluster- (level 2) covariate. The degrees of freedom adjustment
therefore uses `number_controls()`, which counts all covariates in
`control_formula`. This is equivalent to counting level-2 covariates
only, matching the technical supplement's benchmark degrees of freedom
of `J - K - g - 1` (`g` = number of level-2 covariates). Contrast with
[`design_based_estimators_individual()`](https://lmiratrix.github.io/clusterRCT/reference/design_based_estimators_individual.md),
which works on individual-level data where `control_formula` may mix
level-1 and level-2 covariates, and so explicitly uses
`number_level2_controls()` to count only the level-2 ones.
