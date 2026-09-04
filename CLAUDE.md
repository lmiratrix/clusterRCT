# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project status

This is `clusterRCT`, an R package (companion codebase for a research
project) comparing many different estimation strategies for the Average
Treatment Effect (ATE) in blocked, cluster-randomized trials (CRTs). It
is under active development and is not a polished CRAN package —
README.md explicitly warns it is not stable for outside use. The package
depends on a second GitHub package, `blkvar`
(`devtools::install_github("lmiratrix/blkvar")`), which is not on CRAN.

## Common commands

This is a standard `devtools`/`testthat` (edition 3) R package. Run
these from an R session with the working directory set to the package
root (or open `clusterRCT.Rproj` in RStudio).

- Load package for interactive dev: `devtools::load_all()`
- Run full test suite: `devtools::test()`
- Run a single test file:
  `devtools::test_file("tests/testthat/test-design_based_estimators.R")`
- Run tests matching a name/description:
  `devtools::test(filter = "aggregation")`
- Regenerate docs/NAMESPACE from roxygen comments after editing
  `@export`/`@param` etc.: `devtools::document()`
- Full check (build + tests + docs): `devtools::check()`
- Build the vignette: `devtools::build_vignettes()` (source at
  `vignettes/overview_estimators.Rmd`)

There is no separate linter config; style is enforced by convention (see
below), not by a CI lint step.

## Core data model: the canonical formula

Nearly every function in this package accepts data via a formula of the
form:

    Yobs ~ Z | clusterID          # unblocked
    Yobs ~ Z | clusterID | blockID  # blocked (randomization blocks nested in blockID)

`Z` must be binary (0/1) treatment assignment; `clusterID` is the
randomized unit (e.g., school); `blockID` is an optional stratum.
`deconstruct_var_formula()` (R/helper_functions.R) parses this.
`make_canonical_data()` converts *any* input (formula+data, or a formula
whose LHS is already a data frame) into a canonical internal data frame
with fixed column names `Yobs`, `Z`, `clusterID`, `blockID` (and control
covariates appended as-is) — nearly all internal estimator functions
take `formula = NULL` and operate directly on canonical data rather than
re-parsing a user formula. When adding a new estimator function, follow
this same pattern: accept `formula`/`data`/`control_formula`, call
`make_canonical_data()` once at the top, then operate on canonical
columns.

Control covariates are passed separately via `control_formula`
(e.g. `~ X1 + X2`, no LHS); `deconstruct_control_formula()` validates
these don’t collide with the reserved names `Y`, `Z`, `clusterID`,
`blockID`.

## Estimator families and naming convention

[`compare_methods()`](https://lmiratrix.github.io/clusterRCT/reference/compare_methods.md)
(R/compare_methods.R) is the main entry point: given a formula/data, it
calls out to each estimator family, row-binds the results, and returns
one data frame with a point estimate + SE per method. Each family lives
in its own file and is independently callable:

- `R/linear_model_estimators.R` — OLS-based estimators: **LR**
  (person-weighted linear regression), **AR**
  (aggregated/cluster-weighted regression);
  `interacted_linear_model_estimators` adds block-by-treatment
  interactions for blocked designs.
- `R/helper_interacted_estimators.R` — shared machinery for combining
  block-level interacted estimates into an overall ATE
  (`generate_all_interacted_estimates`, precision weighting via
  `calc_agg_estimate`, Welch-Satterthwaite df).
- `R/aggregation_estimators.R` — estimators run on cluster-aggregated
  data (see
  [`aggregate_data()`](https://lmiratrix.github.io/clusterRCT/reference/aggregate_data.md)).
- `R/design_based_estimators.R` — **DB_HT** / **DB_Raj** design-based
  estimators following Schochet, Pashley, Miratrix & Kautz (JASA), plus
  the Middleton & Aronow Raj-type estimator.
- `R/MDStdCRT_estimator.R` — **MRStdCRT** standardization-based
  estimator (see also the standalone `MRStdCRT/` and
  `trial_of_MRStdCRT.R` exploratory area at repo root, not yet
  integrated into the package proper).
- `R/MLM_model_estimators.R` — multilevel/mixed-effects model (**MLM**)
  estimators via `lme4`/`lmerTest`, with fixed-effect (**FE**),
  random-intercept (**RE**), and fixed/random-interaction
  (**FIRC**/**RIRC**) block-handling variants.
- `R/gee_estimators.R` — **GEE** estimators via `geepack`.
- `R/canonical_weight_esimators.R` — canonical-weighting estimators
  (note the filename typo — it’s `esimators.R`, not `estimators.R`).

Method names encode their configuration as dash-separated suffixes,
e.g. `LR-FIcw-db` = linear regression,
fixed-intercept-by-cluster-weighted block handling, design-based SEs.
[`method_characteristics()`](https://lmiratrix.github.io/clusterRCT/reference/method_characteristics.md)
is the canonical lookup table of every implemented method name to its
properties (regression weighting scheme — Person vs. Cluster — whether
it’s biased, whether it requires blocking, and whether it’s “disfavored”
due to odd weighting/instability).
[`get_estimand()`](https://lmiratrix.github.io/clusterRCT/reference/get_estimand.md)
looks up the target estimand (population weighting) for a given method
name. When adding a new method variant, register its
name/characteristics in
[`method_characteristics()`](https://lmiratrix.github.io/clusterRCT/reference/method_characteristics.md)
(R/compare_methods.R) so it flows through
[`get_estimand()`](https://lmiratrix.github.io/clusterRCT/reference/get_estimand.md)
and any reporting that filters on `include_disfavored`.

## Other key pieces

- `R/describe_clusterRCT.R` —
  [`describe_clusterRCT()`](https://lmiratrix.github.io/clusterRCT/reference/describe_clusterRCT.md)
  summarizes structural characteristics of a dataset (block/cluster
  sizes, ICCs via `calc_ICCs()`, nesting via
  [`is_nested()`](https://lmiratrix.github.io/clusterRCT/reference/is_nested.md));
  returns a `clusterRCTstats` S3 object with `print`/`as.data.frame`
  methods.
- `R/patch_data_set.R` —
  [`patch_data_set()`](https://lmiratrix.github.io/clusterRCT/reference/patch_data_set.md)
  handles missing data (mean-imputation with dummy flags) and
  [`patch_singleton_blocks()`](https://lmiratrix.github.io/clusterRCT/reference/patch_singleton_blocks.md)
  handles blocks/clusters with no variation in treatment (all-treated or
  all-control), which break variance estimation;
  [`compare_methods()`](https://lmiratrix.github.io/clusterRCT/reference/compare_methods.md)’s
  `handle_singleton_blocks` argument (`"drop"`/`"pool"`/`"fail"`)
  controls this.
- `R/rerandomize_simulation_code.R` —
  [`rerandomize()`](https://lmiratrix.github.io/clusterRCT/reference/rerandomize.md)
  re-simulates the randomization (respecting cluster/block structure)
  for randomization-inference-style simulation studies;
  [`run_rerandomize_simulation()`](https://lmiratrix.github.io/clusterRCT/reference/run_rerandomize_simulation.md)
  /
  [`summarize_simulation_results()`](https://lmiratrix.github.io/clusterRCT/reference/summarize_simulation_results.md)
  drive Monte Carlo comparisons of the estimators above under known
  ground truth.
- Test data: `fakeCRT`, `fakeCRT2`, `fakeBrokeCRT` (in `data/`,
  generated via `data-raw/fakeCRT.R` using the PUMP package) are 3-level
  (student-in-school-in-district) simulated datasets used throughout the
  tests and vignette; `fakeBrokeCRT` intentionally has missing data and
  singleton blocks to exercise the patching logic.

## Non-package directories at repo root

`MRStdCRT/`, `clusterRCT drafts etc/`, and `luke_scratch/` are
exploratory/scratch work not part of the installable package (the latter
is gitignored). Don’t assume code there is exported, tested, or
maintained to the same standard as `R/`.
