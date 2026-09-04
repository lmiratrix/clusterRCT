# Describe a cluster RCT dataset

Calculate various statistics describing total clusters, etc. If formula
is NULL, data assumed to be in canonical form.

## Usage

``` r
describe_clusterRCT(
  formula = NULL,
  data = NULL,
  control_formula = NULL,
  warn_missing = TRUE
)
```

## Arguments

- formula:

  Formula of the form \`Yobs ~ Z \| clusterID\` (or \`Yobs ~ Z \|
  clusterID \| blockID\` if blocked). If NULL, \`data\` is assumed to
  already be in canonical form (columns \`Yobs\`, \`Z\`, \`clusterID\`,
  and optionally \`blockID\`).

- data:

  Data frame holding the outcome, treatment, and clustering/blocking
  columns referenced by \`formula\` (or already in canonical form if
  \`formula\` is NULL).

- control_formula:

  What variables to control for, in the form of "~ X1 + X2".

- warn_missing:

  If TRUE, warn when rows are dropped or imputed due to missing data.

## Value

A list with the following elements:

- n:

  Total number of individuals.

- J:

  Total number of clusters.

- nbar:

  Average cluster size (n/J).

- ncv:

  Coefficient of variation of cluster sizes (sd(n_j)/nbar).

- n.25:

  25th percentile of cluster sizes.

- n.75:

  75th percentile of cluster sizes.

- n.IQR:

  Interquartile range of cluster sizes.

- p.tx:

  Proportion of clusters assigned to treatment.

- K:

  Total number of blocks.

- Jbar:

  Average number of clusters per block.

- Jcv:

  Coefficient of variation of clusters per block.

- J.25:

  25th percentile of clusters per block.

- J.75:

  75th percentile of clusters per block.

- J.IQR:

  Interquartile range of clusters per block.

- n_block:

  Average number of individuals per block.

- n_block_cv:

  Coefficient of variation of individuals per block.

- tx.avg:

  Average proportion of clusters treated across blocks (not the same as
  p.tx).

- tx.cv:

  Coefficient of variation of proportion clusters treated across blocks.

- tx.25:

  25th percentile of proportion clusters treated across blocks.

- tx.75:

  75th percentile of proportion clusters treated across blocks.

- tx.IQR:

  Interquartile range of proportion clusters treated across blocks.

- cluster_ICC:

  Intraclass correlation coefficient at the cluster level.

- block_ICC:

  Intraclass correlation coefficient at the block level.

- sdY0:

  Estimated standard deviation of the outcome under control.

- matched_pairs:

  TRUE/FALSE of being matched pairs experiment.

- num_singletons:

  Number of treatment arms wihtin a block with only a single cluster.

- num_doubletons:

  Number of treatment arms wihtin a block with exactly two clusters.

## Details

Will patch data via patch_data() to handle missing values.

## See also

\[patch_data_set()\]
