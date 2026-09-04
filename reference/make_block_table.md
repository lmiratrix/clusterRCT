# Make table of statistics by block

Given individual/school/block data of a blocked, cluster RCT, calculate
statistics for each block (block).

## Usage

``` r
make_block_table(
  formula = NULL,
  data = NULL,
  control_formula = NULL,
  check_data_integrity = FALSE
)
```

## Arguments

- formula:

  Notation for Y ~ Z \| clusterID \| blockID. If NULL, \`data\` is
  assumed to already be in canonical form.

- data:

  Dataframe holding the outcome, treatment, and clustering/blocking
  columns referenced by \`formula\` (or already in canonical form if
  \`formula\` is NULL).

- control_formula:

  What variables to control for, in the form of "~ X1 + X2".

- check_data_integrity:

  TRUE means runs some checks and give errors if data fails them (e.g.,
  incorrectly processed treatment vector.). FALSE means calculate
  statistics without these checks.
