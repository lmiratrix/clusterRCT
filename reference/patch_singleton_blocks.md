# Drop or combine all-tx or all-co blocks

Given dataset with some blocks that are all treated or all control, drop
those blocks or replace the block ID with canonical new block ID shared
by all such blocks.

## Usage

``` r
patch_singleton_blocks(
  formula = NULL,
  data,
  drop_data = TRUE,
  pool_clusters = TRUE,
  warn_missing = TRUE
)
```

## Arguments

- formula:

  Notation for Y ~ Z \| clusterID \| blockID. If NULL, \`data\` is
  assumed to already be in canonical form.

- data:

  Dataframe to patch.

- drop_data:

  Drop the troublesome blocks if TRUE, pool them if FALSE.

- pool_clusters:

  If pooling blocks rather than dropping them, also pool clusters in 100
  cluster.

- warn_missing:

  Say something if anything happens.

## Value

data with the blockID column modified or troublesome rows of data
missing.

## Details

Missing values in the block ID are all considered singletons.

Also, depending on pool_clusters, pool the clusters in each of these
identified blocks into single clusters.
