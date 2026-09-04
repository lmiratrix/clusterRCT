# Run a finite sample simulation by permuting treatment assignment labels within each block.

Before starting simulation, put data in canonical form, imputing missing
values and so forth as needed. Call compare_methods repeatidly on the
resulting treatment-permuted data.

## Usage

``` r
run_rerandomize_simulation(
  formula,
  data,
  R = 100,
  control_formula = NULL,
  summarize_results = FALSE,
  include_empirical = TRUE,
  warn_missing = TRUE,
  patch_data = TRUE,
  parallel = FALSE,
  ...
)
```

## Arguments

- formula:

  A formula object specifying the model to estimate.

- data:

  A data frame containing the data to analyze.

- R:

  The number of simulations to run.

- control_formula:

  A formula object specifying the control variables; passed to
  compare_methods().

- summarize_results:

  If TRUE, summarize the results of the compare_methods call, FALSE
  return raw estimates.

- include_empirical:

  If TRUE, include the empirical results as one of the "simulation
  replicates". This will be included in any summarization, so be warned.

- warn_missing:

  If TRUE, warn when rows are dropped or covariates are imputed while
  canonicalizing/patching the data.

- patch_data:

  If TRUE, impute missing covariates via \`patch_data_set()\` before
  simulating.

- parallel:

  If TRUE, run simulations in parallel.

- ...:

  Additional arguments to pass to compare_methods.
