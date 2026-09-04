# Cluster RCT with multi-level modeling

This uses the lmerTest package to get the p-values, and fits a model
with assumed homoskedasticity, etc. I.e., this is the vanilla MLM that
one would typically fit.

## Usage

``` r
MLM_estimators(
  formula,
  data = NULL,
  control_formula = NULL,
  suppress_warnings = TRUE,
  include_disfavored = FALSE
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

- suppress_warnings:

  If TRUE, suppress the (typically benign) convergence/singularity
  warnings \`lmer()\` throws.

- include_disfavored:

  Include the full random-intercept- and random-slope-by-block
  (RIRC/FIRC) MLM variants, which are considered disfavored due to
  instability with few blocks.
