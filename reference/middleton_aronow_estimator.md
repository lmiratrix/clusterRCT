# Estimate ATE via Middleton & Aronow unbiased approach

This follows Middleton, J. A., & Aronow, P. M. (2015). Unbiased
Estimation of the Average Treatment Effect in Cluster-Randomized
Experiments. 6(1–2), 312. doi: 10.1515/spp-2013-0002

## Usage

``` r
middleton_aronow_estimator(
  formula,
  data = NULL,
  control_formula = NULL,
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

They introduce a Raj Difference estimator that is an extension of a
Horvitz-Thompson estimator where, in effect, we divide the sum of the
outcomes by a fixed constant rather than the realized sample size to
avoid biasing our estimate.

Note that the control_formula is ignored–these estimators do not adjust
for additioanl controls at this time.
