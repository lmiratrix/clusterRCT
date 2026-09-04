# Make a canonical regression formula, possibly with control variables.

Make a canonical regression formula, possibly with control variables.

## Usage

``` r
make_regression_formula(
  Yobs = "Yobs",
  Z = "Z",
  clusterID = "clusterID",
  blockID = "blockID",
  control_formula = NULL,
  interacted = FALSE,
  control_interacted = FALSE,
  FE = FALSE,
  cluster_RE = FALSE,
  data = NULL
)
```

## Arguments

- Yobs:

  Name of outcome variable (assumed to exist in data)

- Z:

  Name of treatment variable (assumed to exist in data)

- clusterID:

  Name of cluster ID variable (assumed to exist in data)

- blockID:

  Name of block ID variable (assumed to exist in data, if this is not
  null).

- control_formula:

  What variables to control for, in the form of "~ X1 + X2".

- interacted:

  TRUE means include treatment by block interactions. Will override FE
  flag and set FE to true.

- control_interacted:

  TRUE means include treatment by control interactions in regression

- FE:

  TRUE means include block dummy variables.

- cluster_RE:

  Add a random effect term for clusterID

- data:

  Dataframe holding all variables to be used in formula.

## Value

Something like "Yobs ~ 1 + Z" or "Yobs ~ 1 + Z + X1 + X2"
