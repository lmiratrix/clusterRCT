# clusterRCT public-release checklist

Generated from a full code review (2026-09-03/04), worked through after the open questions were answered (2026-09-04). Work happened on a git branch, `cleanup/public-release-prep`, off `main` -- **main is untouched**. The commit that produced the accepted JREE paper's simulation results is tagged `jree-submission-snapshot-2025-11-20` (confirmed by you).

Deferred items needing statistical verification live in `FUTURE_WORK.md`, not here.

**Overall result: `devtools::check()` (including a full vignette rebuild) now passes with 0 errors, 0 warnings, 1 note** (down from 7 warnings + 4 notes at the start). Full `devtools::test()` suite passes throughout with 0 failures.

## Packaging / dependency fixes

- [x] Added missing runtime dependencies to DESCRIPTION: `geepack`, `broom`, `stringr` (all used but previously undeclared -- would crash for any user without them already installed).
- [x] Added `furrr`/`future`/`tidyverse`/`PUMP` to Suggests (optional, dev-only, or test-only).
- [x] `run_rerandomize_simulation(parallel = TRUE)` now fails with a clear "install furrr" message instead of a bare `library(furrr)` crash.
- [x] Removed 10 declared-but-unused Imports (`Formula`, `MASS`, `ggplot2`, `lmtest`, `magrittr`, `msm`, `nlme`, `rlang`, `sandwich`, `survey`).
- [x] Added `@importFrom lme4 VarCorr fixef lmer` and `@importFrom tidyr pivot_wider` (fixes "Package in Depends field not imported from" and removes the need for a `require(tidyr)` workaround).
- [x] Added `@importFrom stats ...` / `@importFrom utils hasName` for ~18 base functions used unqualified (`as.formula`, `coef`, `lm`, `model.matrix`, `sd`, `quantile`, `pnorm`, `pt`, `var`, `vcov`, `weighted.mean`, etc.).
- [x] Removed dead/redundant in-body `require()`/`library()` calls that were 100% no-ops: `require(dplyr)`, `require(lme4)`, `require(lmerTest)` (x2 files), `require(formula.tools)` (x3 sites), `require(estimatr)` (x5 sites), `require(MRStdCRT)`.
- [x] Replaced unnecessary self-referential `clusterRCT:::foo()` calls with direct calls (5 files).
- [x] Rewrote `LICENSE` as the proper short DCF stub; full text moved to `LICENSE.md`.
- [x] Added `MRStdCRT/`, `clusterRCT drafts etc/`, `luke_scratch/`, `rctYEScheck/`, `trial_of_MRStdCRT.R`, `Mike's Notes on Output.docx`, `.claude` to `.Rbuildignore` (this is what actually controls the built package -- `.gitignore` does not) and to `.gitignore`.

## Bugs fixed

- [x] **`calc_ICCs()` unblocked-branch normalization bug**: `C.ICC` was left unnormalized (typo'd into a dead variable `C.ISS` instead of updating `C.ICC`) instead of `C.ICC / tvar`. Fixed; verified `describe_clusterRCT()` on unblocked data now returns an ICC in \[0,1\].
- [x] `as.data.frame.clusterRCTstats()`'s S3 method signature (was missing `row.names`/`optional`/`...`).
- [x] `compare_methods()`'s multi-outcome recursion was silently dropping `include_gee`/`include_dumb`/`include_disfavored` on recursive per-outcome calls -- now passes through correctly.
- [x] Duplicate Rd alias (`clusterRCTstats` declared on two different `.Rd` files).
- [x] Stale/incomplete data documentation: `fakeCRT`/`fakeBrokeCRT`'s undocumented `X` column, `fakeCRT2`'s undocumented `X.jk`/`C.ijk` columns and wrong row count (860, not 1500).
- [x] `canonical_weight_estimators()`'s `VarCorr(...)[i,j]` positional indexing replaced with name-based indexing -- verified against actual `lme4` output structure before changing.
- [x] `method_characteristics()`/`compare_methods()` silent-drop: methods not registered in `method_characteristics()` used to vanish entirely from output. Now kept with `NA` characteristics.

## API changes (confirmed with you)

- [x] `gee_estimators()`: removed `"Cluster"` from the documented/selectable `weight` options; left the `stop()` + unreachable weighting skeleton in the code for future development.
- [x] `compare_methods()`: default `include_gee` changed from `TRUE` to `FALSE` (GEE is still developing and can hit rank-deficient designs -- this is also what was silently breaking the vignette build).
- [x] `compare_methods()`: removed the unused `weight` parameter entirely.
- [x] `MRStdCRT_estimator()`: unexported (`@export` -> `@noRd`), `require(MRStdCRT)` removed. Left non-functional as a skeleton for future work -- see `FUTURE_WORK.md`.
- [x] `canonical_weight_estimators()`: unexported (`@export` -> `@noRd`) since it's a diagnostic/validation-only tool, not a real estimator for end users.

## Documentation

- [x] `design_based_estimators()` / `design_based_estimators_individual()`: added `@details` explaining which covariate-count helper each uses for its degrees-of-freedom adjustment and why they're expected to differ. See `FUTURE_WORK.md` for the "confirm this is exactly right" follow-up.
- [x] Full `@param` documentation sweep: every exported function now passes `devtools::check()`'s "Rd \\usage sections" check (was ~24 functions with undocumented parameters, including stale `Yobs`/`Z`/`B`/`blockID` params left over from an old pre-formula `compare_methods()` signature). Several internal-only helper functions (`calc_ICCs`, `calc_covariate_R2s`, `check_data_integrity`, `deconstruct_control_formula`, `deconstruct_var_formula`, `expand_control_variables`) were `@noRd`'d instead, since they were never exported and don't need public-facing docs.
- [x] Removed the "USERS TAKE NOTE ... I would not use it if I were you" disclaimer from README.md.
- [x] Removed the stale `blkvar` GitHub install instruction from README.md (confirmed no longer needed); replaced with a real `devtools::install_github("lmiratrix/clusterRCT")` instruction.
- [x] Fixed the vignette's leftover "this doesn't work yet since package is private" commented-out install block; vignette now builds cleanly end-to-end.

## `.data`/NSE cleanup ("no visible binding for global variable" NOTE)

Went from ~200 flagged bindings across ~30 functions down to ~20, across every live (non-`if(FALSE)`) code path in `compare_methods.R`, `describe_clusterRCT.R`, `helper_functions.R`, `patch_data_set.R`, `design_based_estimators.R` (including the dense Schochet variance-formula functions -- converted carefully and verified against the file's hand-computed-expected-value test), `linear_model_estimators.R`, `aggregation_estimators.R`, `gee_estimators.R`, `canonical_weight_esimators.R`, and `rerandomize_simulation_code.R`.

- [x] Converted bare column references inside `mutate()`/`summarise()`/`filter()`/`group_by()`/`arrange()` to the `.data$col` pronoun (verified empirically that this produces byte-identical output to bare references).
- [x] Converted the equivalent references inside tidyselect verbs (`select()`, `relocate()`, `rename()`) to string literals instead, after discovering `.data$col` is deprecated specifically in tidyselect contexts (caught via a live deprecation warning during testing, fixed, and reverified).
- [x] Where `data$col` couldn't be substituted (base `stats::lm()`'s `weights=`/`estimatr::lm_robust()`'s NSE evaluates certain arguments in a way that broke across function-call frames -- caught by an actual test failure, then correctly reverted and diagnosed), left the original bare-symbol form in place rather than risk a subtle correctness bug in the paper's core SE formulas. See below.
- [x] Every single change in this section was verified against the full test suite (0 failures throughout, including a live catch-and-fix of a base-`lm()` NSE frame issue in `design_based_estimators.R` that would have silently broken standard error calculations had it not been caught by `test-design_based_estimators.R`'s hand-computed-expected-value tests).

**Remaining ~20 bindings, deliberately left as-is** (not bugs, just NOTE-level cosmetics that don't have a safe `.data$` fix):
- `MRStdCRT_estimator`'s dead/unreachable code (future-work skeleton, not touched).
- `lm(weights = .weight)` (6 call sites in `design_based_estimators.R`) and `geeglm(id = clusterID)` (1 site in `gee_estimators.R`): base-R-style modeling functions whose NSE evaluates `weights=`/`clusters=`/`id=` arguments against the `data=` argument's own columns, not the calling frame -- passing `data$col` directly breaks this (confirmed by a real test failure) for `stats::lm()`, so the bare form is intentionally kept. This is a well-known, common, and unavoidable NOTE for any package using `lm(weights=...)`.
- A handful of ambiguous cases (`m`, `method`, `weight` in `generate_all_interacted_estimates()`/`canonical_weight_estimators()`/`compare_methods()`) where the bare name is a local function parameter or scalar variable, not a data column -- `.data$` would be actively wrong here (would error or silently look up the wrong thing), so left alone.

## Test coverage added

- [x] `tests/testthat/test-rerandomize.R`: `rerandomize()` had zero direct test coverage despite being the foundation of every randomization-inference simulation. Added tests for cluster-constant treatment, block/overall proportion preservation, and a Monte Carlo check.
- [x] `tests/testthat/test-canonical_weight_estimators.R`: added a real correctness test verifying the 6 non-MLM canonical weight combinations exactly reproduce their corresponding linear-model estimator's point estimate.
- [x] Fixed a latent test-isolation bug in `test-patch_data_set.R` (relied on an unqualified `rhs.vars()` call that only worked by accident via a since-removed `require(formula.tools)` side effect elsewhere in the package).

## Still open (see `FUTURE_WORK.md` for the statistical items)

- [ ] No CI (`.github/workflows`) -- worth an `R-CMD-check` GitHub Action now that `devtools::check()` is clean.
- [ ] No `NEWS.md`, no `URL`/`BugReports` fields in DESCRIPTION, no `cran-comments.md`, no pkgdown site.
- [ ] Test coverage still missing for: `design_based_estimators_individual()` (zero), `MLM_estimators(include_disfavored = TRUE)` RIRC/FIRC path (zero), `aggregation_estimators()` (smoke-tested only), GEE's unblocked/`control_formula` paths (untested).
- [ ] Personal/dev directories at repo root (`MRStdCRT/`, `clusterRCT drafts etc/`, `luke_scratch/`, `rctYEScheck/`, `trial_of_MRStdCRT.R`) are now excluded from the built package via `.Rbuildignore` but still sitting in the repo -- decide whether to delete/relocate/leave.
- [ ] `import(dplyr)` + `Depends: dplyr` is still a blanket import rather than selective `importFrom` -- not a functional bug (verified), just not best practice; low priority.
