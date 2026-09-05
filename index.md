# clusterRCT - code for analyzing blocked, cluster-randomized trials

# Introduction

`clusterRCT` provides an easy means to run a collection of different
estimators to estimate average treatment effects in blocked, cluster
randomized trials (trials where you have randomization blocks, and
within each block clusters are randomized into treatment and control).
It also provides some methods for exploring data with a blocked,
clustered structure. The package also can be used for simple
cluster-randomized trials with no blocking.

Core methods are
[`compare_methods()`](https://lmiratrix.github.io/clusterRCT/reference/compare_methods.md)
which compares the estimates from a collection of different estimation
strategies on a given dataset, and
[`describe_clusterRCT()`](https://lmiratrix.github.io/clusterRCT/reference/describe_clusterRCT.md)
which describes overall characteristics of a given dataset.

The package was designed to support a paper, “An Applied Researcher’s
Guide to Estimating Effects from Multisite Cluster Randomized Trials:
Estimands, Estimators, and Estimates,” accepted to JREE in 2026. For the
replication files of the paper, see [the replication
GitHub](https://github.com/lmiratrix/clusterRCTsim). That set of code
heavily uses this package in the implemented simulations and empirical
data analyses.

# Installation

Installation is simple:

    devtools::install_github("lmiratrix/clusterRCT")

## Replicating the paper

The paper’s results were generated with a slightly older, less
cleaned-up version of the package than what is on the main branch today.
Ongoing cleanup (dependency fixes, bug fixes, documentation) is not
intended to change any reported estimates, but to install the exact
version used to produce the paper’s results, use the tagged release
instead of the latest version:

    devtools::install_github("lmiratrix/clusterRCT",
                ref = "jree-submission-snapshot-2025-11-20")

## Acknowledgements

The research reported here was partially supported by the Institute of
Education Sciences, U.S. Department of Education, through Grant
R305D220046. The opinions expressed are those of the authors and do not
represent views of the Institute or the U.S. Department of Education.
