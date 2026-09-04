# Estimate ATEs for a cluster RCT using a GEE approach

Estimate ATEs for a cluster RCT using a GEE approach

## Usage

``` r
gee_estimators(
  formula,
  data = NULL,
  control_formula = NULL,
  weight = c("Person")
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

  Weighting scheme for the GEE. Only "Person" is currently implemented;
  see FUTURE_WORK.md for cluster weighting.
