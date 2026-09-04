# Fake Broken Cluster Randomized Trial data

This dataset is a version of \`fakeCRT\` with missing data and some
blocks that are all tx or all co.

## Usage

``` r
fakeBrokeCRT
```

## Format

A data frame with 1500 rows and 8 variables:

- `V.k`:

  double District-level covariate (think standardized district-wide
  average test score relative to state).

- `X.jk`:

  double School-level (think percent on Free/Reduced Price Lunch).

- `C.ijk`:

  double Student-level covariate (think baseline measured SES).

- `S.id`:

  integer School ID.

- `D.id`:

  integer District ID.

- `Yobs`:

  double Observed outcome (think test score).)

- `T.x`:

  integer Treatment assignment (1 treated, 0 control).

- `X`:

  character A categorical (4-level) student-level covariate, for testing
  categorical covariate handling.

## Details

These data were generated via the PUMP package.
