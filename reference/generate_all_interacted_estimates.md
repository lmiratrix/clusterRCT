# Average block-level estimates to get overall ATE

Also calculate associated standard error.

## Usage

``` r
generate_all_interacted_estimates(
  fitModel,
  data,
  use_full_vcov = FALSE,
  SE_table = NULL,
  method = "LR",
  weight = NULL,
  se_method = "crve",
  aggregated = FALSE,
  include_block_estimates = FALSE
)
```

## Arguments

- fitModel:

  The model that has the interacted estimates in it.

- data:

  Data frame the model was fit on (in canonical form).

- use_full_vcov:

  TRUE/FALSE. When calculating standard errors, should the full
  variance-covariance matrix of the block-level estimates be used, or
  just the diagonal of standard errors?

- SE_table:

  Optional precomputed table of block-level standard errors to use
  instead of extracting them from \`fitModel\`.

- method:

  Prefix of the method. Will add the weighting to stem.

- weight:

  "Person" or "Cluster" weighting used when averaging block-level
  estimates to an overall ATE.

- se_method:

  Which standard error method the block-level estimates in
  \`fitModel\`/\`SE_table\` were computed with (e.g. "crve", "het",
  "db"). Used only for labeling.

- aggregated:

  TRUE means data is summarized at the Block level. False means it is
  not, and needs to be aggregated.

- include_block_estimates:

  If TRUE, also return the block-by-block estimates used to build the
  overall estimate, in addition to the aggregate estimate.
