# Check if one grouping variable nested in another

Checks if clusterID nested in blockID, so each clusterID value appears
in only one blockID value.

## Usage

``` r
is_nested(clusterID, blockID)
```

## Arguments

- clusterID:

  List of clusterIDs

- blockID:

  List of blockIDs

## Value

TRUE if nested, FALSE otherwise. If blockID NULL, always return TRUE.
