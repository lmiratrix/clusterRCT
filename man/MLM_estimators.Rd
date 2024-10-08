% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLM_model_estimators.R
\name{MLM_estimators}
\alias{MLM_estimators}
\title{Cluster RCT with multi-level modeling}
\usage{
MLM_estimators(
  formula,
  data = NULL,
  control_formula = NULL,
  suppress_warnings = TRUE
)
}
\arguments{
\item{formula}{Formula for outcome and treatment and nesting.  If
NULL, data is assumed to be in canonical form (see vignette for
further discussion).}
}
\description{
This uses the lmerTest package to get the p-values, and fits a
model with assumed homoskedasticity, etc.  I.e., this is the
vanilla MLM that one would typically fit.
}
