# clusterRCT - code for analyzing blocked, cluster-randomized trials

# Introduction
This is the companion codebase for a large-scale comparison of different estimators applied to a suite of cluster RCT experiments in Education and the social sciences.

Core methods are `compare_methods()` which compares the estimates from a collection of different estimation strategies on a given dataset, and `describe_clusterRCT()` which describes overall characteristics of a given dataset.


# Installation

```
devtools::install_github("lmiratrix/clusterRCT")
```

## Replicating the paper

This package is the companion to "An Applied Researcher's Guide to Estimating Effects from Multisite Cluster Randomized Trials: Estimands, Estimators, and Estimates." For the paper's replication files, we used a slightly older, less cleaned-up version of the package than what is on the main branch today. Ongoing cleanup (dependency fixes, bug fixes, documentation) is not intended to change any reported estimates, but to install the exact version used to produce the paper's results, use the tagged release instead of the latest version:

```
devtools::install_github("lmiratrix/clusterRCT", ref = "jree-submission-snapshot-2025-11-20")
```


## Acknowledgements
The research reported here was partially supported by the Institute of Education Sciences, U.S. Department of Education, through Grant R305D220046.
The opinions expressed are those of the authors and do not represent views of the Institute or the U.S. Department of Education.

## Fake Version Tag

The version of this code is at least 0.006 as of March 26, 2024



