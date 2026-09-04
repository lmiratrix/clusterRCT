# Identify all tx and all co blocks

Identify blocks that are all tx or all co. Rows with missing treatment
indicators are dropped. Rows with missing block identifiers are
considered all unique blocks and thus put into the table of results.

## Usage

``` r
identify_singleton_blocks(formula, data)
```

## Arguments

- formula:

  Notation for Y ~ Z \| clusterID \| blockID. If NULL, \`data\` is
  assumed to already be in canonical form.

- data:

  Dataframe to check.

## Value

Dataframe of block IDs corresponding to all tx or all co blocks along
with number of units and tx status. Returns NULL if nothing.
