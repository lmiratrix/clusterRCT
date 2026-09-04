# Estimating ATEs for CRTs

This package is designed to make comparing different estimators for the
average treatment effect (ATE) of a cluster randomized trial (CRT)
relatively straightforward. There are a host of methods that will
estimate the ATE using different modeling approaches (linear regression
models, multilevel models, design-based methods, aggregate-then-analyze
approaches). There is also an overall method
[`compare_methods()`](https://lmiratrix.github.io/clusterRCT/reference/compare_methods.md)
that compares all the estimators in one large table of results.

To install the package:

``` r

devtools::install_github( "lmiratrix/clusterRCT" )
```

Once installed, to load the package, first call
[`library()`](https://rdrr.io/r/base/library.html):

``` r

library(clusterRCT)
```

To illustrate we use a fake dataset embedded in the package called
`fakeCRT`:

``` r

data( fakeCRT )
head( fakeCRT )
#> # A tibble: 6 × 8
#>       V.k   X.jk   C.ijk S.id  D.id    Yobs   T.x X    
#>     <dbl>  <dbl>   <dbl> <fct> <fct>  <dbl> <int> <chr>
#> 1  1.53   -0.988 -0.0353 164   6      0.499     0 A    
#> 2 -1.69   -1.10  -0.161  192   7     -0.769     0 A    
#> 3 -1.69   -1.10  -1.77   192   7     -1.33      0 B    
#> 4  1.45    0.735 -0.675  134   5     -1.33      0 D    
#> 5  1.53    0.666  1.14   175   6     -0.979     0 B    
#> 6 -0.0386 -0.485 -0.628  80    3      0.244     0 A
```

We have a few covariates, a treatment assignment, and an outcome. The
`S.id` and `D.id` are school and district IDs, respectively.

We can estimate impacts as follows:

``` r

compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT, 
                 include_method_characteristics = FALSE ) %>%
    knitr::kable( digits = 2 )
```

| method         | ATE_hat | SE_hat | p_value |     df |
|:---------------|--------:|-------:|--------:|-------:|
| AR-FE-db       |    0.19 |   0.09 |    0.04 | 189.00 |
| AR-FE-het      |    0.19 |   0.09 |    0.04 | 189.00 |
| AR-FIcw-db     |    0.19 |   0.09 |    0.03 | 180.00 |
| AR-FIcw-het    |    0.19 |   0.09 |    0.03 | 180.00 |
| ARpw-FE-db     |    0.20 |   0.09 |    0.02 | 189.00 |
| ARpw-FE-het    |    0.20 |   0.09 |    0.02 | 189.00 |
| ARpw-FIpw-db   |    0.20 |   0.09 |    0.02 | 180.00 |
| ARpw-FIpw-het  |    0.20 |   0.09 |    0.02 | 180.00 |
| LR-FE-crve     |    0.20 |   0.09 |    0.02 | 157.66 |
| LR-FE-db       |    0.20 |   0.09 |    0.02 | 189.00 |
| LR-FIpw-crve   |    0.20 |   0.09 |    0.02 | 140.22 |
| LR-FIpw-db     |    0.20 |   0.09 |    0.02 | 180.00 |
| LRcw-FE-crve   |    0.19 |   0.09 |    0.04 | 144.95 |
| LRcw-FE-db     |    0.19 |   0.09 |    0.04 | 189.00 |
| LRcw-FIcw-crve |    0.19 |   0.09 |    0.03 | 134.78 |
| LRcw-FIcw-db   |    0.19 |   0.09 |    0.03 | 180.00 |
| MLM-FE         |    0.19 |   0.09 |    0.03 | 186.22 |
| MLM-FIcw       |    0.19 |   0.09 |    0.03 | 180.00 |
| MLM-RE         |    0.18 |   0.09 |    0.04 | 187.42 |

Each row of the output represents a different estimator (defined as a
point estimator and standard error estimator pair). The first column is
our estimate ATE, the second the estimated standard error. We also have
a calculated p-value and some further information about the estimator
itself. The `include_method_characteristics` flag will eventually allow
us to have notes on each method, such as whether it is targeting a
finite or superpopulation estimand, and so forth.

If we want to control for covariates, we can as so:

``` r

compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                 control_formula = ~ V.k + X.jk + C.ijk,
                 include_method_characteristics = FALSE )
#>            method   ATE_hat     SE_hat    p_value       df
#> 1        AR-FE-db 0.1564187 0.07740607 0.04473906 186.0000
#> 2       AR-FE-het 0.1564187 0.07695446 0.04350556 187.0000
#> 3      AR-FIcw-db        NA         NA         NA 177.0000
#> 4     AR-FIcw-het 0.1614133 0.07727308 0.03671960 177.0000
#> 5      ARpw-FE-db 0.1636241 0.07841208 0.03827674 186.0000
#> 6     ARpw-FE-het 0.1636241 0.07718937 0.03534406 187.0000
#> 7    ARpw-FIpw-db        NA         NA         NA 177.0000
#> 8   ARpw-FIpw-het 0.1658707 0.07623035 0.02956166 177.0000
#> 9      LR-FE-crve 0.1621029 0.07753583 0.03821277 152.7281
#> 10       LR-FE-db 0.1621029 0.07841798 0.04009603 187.0000
#> 11   LR-FIpw-crve        NA         NA         NA 141.1301
#> 12     LR-FIpw-db        NA         NA         NA 178.0000
#> 13   LRcw-FE-crve 0.1489494 0.07765916 0.05720770 136.0520
#> 14     LRcw-FE-db 0.1489494 0.07781325 0.05712259 187.0000
#> 15 LRcw-FIcw-crve        NA         NA         NA 135.9229
#> 16   LRcw-FIcw-db        NA         NA         NA 178.0000
#> 17         MLM-FE 0.1539767 0.07667356 0.04608268 184.0841
#> 18       MLM-FIcw        NA         NA         NA 177.0000
#> 19         MLM-RE 0.1480336 0.07656648 0.05471175 185.2041
```

Here we are controlling for district, school, and individual-level
characteristics.

Finally, we can turn on or off different families of estimator. E.g., we
are dropping the aggregation methods in this call:

``` r

compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                 include_agg = FALSE,
                 include_method_characteristics = FALSE )
#>            method   ATE_hat     SE_hat    p_value       df
#> 1      LR-FE-crve 0.2032865 0.08660625 0.02015520 157.6603
#> 2        LR-FE-db 0.2032865 0.08762329 0.02140859 189.0000
#> 3    LR-FIpw-crve 0.2044673 0.08600557 0.01743643 140.2158
#> 4      LR-FIpw-db 0.2044673 0.08567461 0.01700696 180.0000
#> 5    LRcw-FE-crve 0.1852251 0.08710924 0.03516777 144.9548
#> 6      LRcw-FE-db 0.1852251 0.08740229 0.03537707 189.0000
#> 7  LRcw-FIcw-crve 0.1867940 0.08780516 0.03338943 134.7822
#> 8    LRcw-FIcw-db 0.1867940 0.08783710 0.03345371 180.0000
#> 9          MLM-FE 0.1931067 0.08674392 0.02720194 186.2165
#> 10       MLM-FIcw 0.1947340 0.08709019 0.02535151 180.0000
#> 11         MLM-RE 0.1816252 0.08666143 0.03744273 187.4170
```

If we have data with no blocking (district-level grouping) we run like
this:

``` r


compare_methods( Yobs ~ T.x | S.id, data=fakeCRT,
                 control_formula = ~ V.k + X.jk + C.ijk,
                 include_method_characteristics = FALSE )
#>      method    ATE_hat     SE_hat   p_value       df
#> 1     AR-db 0.10707513 0.08799049 0.2251151 195.0000
#> 2    AR-het 0.10707513 0.08903051 0.2305588 195.0000
#> 3   ARpw-db 0.11668887 0.09125813 0.2025336 195.0000
#> 4  ARpw-het 0.11668887 0.09200956 0.2062295 195.0000
#> 5   LR-crve 0.11345422 0.09270457 0.2227040 170.7613
#> 6     LR-db 0.11345422 0.09172813 0.2176215 196.0000
#> 7 LRcw-crve 0.09419665 0.09062817 0.3002902 151.1416
#> 8   LRcw-db 0.09419665 0.08900770 0.2912222 196.0000
#> 9       MLM 0.09729122 0.09069491 0.2847113 195.9409
```
