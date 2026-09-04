# Estimate ATEs for a cluster RCT using linear model

Estimate district-average effects using an interacted FE model and then
aggregate.

## Usage

``` r
interacted_linear_model_estimators(
  formula,
  data = NULL,
  control_formula = NULL,
  weight = c("Person", "Cluster"),
  use_full_vcov = FALSE
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

- use_full_vcov:

  TRUE/FALSE. When calculating standard errors, should the full
  variance-covariance matrix of the block-level estimates be used, or
  just the diagonal of standard errors?
