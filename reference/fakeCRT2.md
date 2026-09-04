# Fake Cluster Randomized Trial data (set 2)

This dataset is a three-level dataset with students in schools in
districts, with schools randomized to treatment and control. It is used
to illustrate and test the clusterRCT pacakge.

## Usage

``` r
fakeCRT2
```

## Format

A data frame with 860 rows and 8 variables:

- `S.id`:

  character School ID.

- `D.id`:

  character District ID.

- `Y1`:

  double Potential outcome under treatment.

- `Y0`:

  double Potential outcome under no treatment.

- `Yobs`:

  double Observed outcome (think test score).)

- `T.x`:

  integer Treatment assignment (1 treated, 0 control).

- `X.jk`:

  double School-level covariate.

- `C.ijk`:

  double Student-level covariate.

## Details

These data were generated via the PUMP package.

## See also

[`fakeCRT`](https://lmiratrix.github.io/clusterRCT/reference/fakeCRT.md)
