# Estimate ATEs for cluster RCT using data aggregation.

This method aggregates data at the cluster level, and analyzes using
linear regression. Can optionally pass pre-aggregated data if desired.

## Usage

``` r
aggregation_estimators(
  formula,
  data = NULL,
  control_formula = NULL,
  control_interacted = FALSE,
  aggregated = FALSE
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

- control_interacted:

  If TRUE, also fit block-by-covariate interacted control models.
  Ignored if \`control_formula\` is NULL.

- aggregated:

  TRUE means \`data\` is already aggregated to the cluster level (and in
  canonical form). FALSE means it is not and will be aggregated
  internally.

## Details

This automatically does both weighting approaches.
