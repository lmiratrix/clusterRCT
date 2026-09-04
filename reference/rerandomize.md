# Return new treatment assignment based on passed one.

Rerandomizes treatment at the cluster level (within block, if
\`blockID\` is given), preserving the original per-block (or overall)
proportion of clusters treated.

## Usage

``` r
rerandomize(Z, clusterID, blockID = NULL)
```

## Arguments

- Z:

  Vector of original treatment assignments (1 = treated, 0 = control),
  one entry per row/unit.

- clusterID:

  Vector of cluster IDs, one entry per row/unit.

- blockID:

  Vector of block IDs, one entry per row/unit. If NULL, rerandomization
  ignores blocking.

## Value

vector of 1s and 0s.
