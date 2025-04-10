% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rerandomize_simulation_code.R
\name{run_rerandomize_simulation}
\alias{run_rerandomize_simulation}
\title{Run a finite sample simulation by permuting treatment assignment
labels within each block.}
\usage{
run_rerandomize_simulation(
  formula,
  data,
  R = 100,
  control_formula = NULL,
  summarize_results = FALSE,
  include_empirical = TRUE,
  warn_missing = TRUE,
  patch_data = TRUE,
  parallel = FALSE,
  ...
)
}
\arguments{
\item{formula}{A formula object specifying the model to estimate.}

\item{data}{A data frame containing the data to analyze.}

\item{R}{The number of simulations to run.}

\item{control_formula}{A formula object specifying the control
variables; passed to compare_methods().}

\item{summarize_results}{If TRUE, summarize the results of the
compare_methods call, FALSE return raw estimates.}

\item{include_empirical}{If TRUE, include the empirical results as
one of the "simulation replicates".  This will be included in any
summarization, so be warned.}

\item{parallel}{If TRUE, run simulations in parallel.}

\item{...}{Additional arguments to pass to compare_methods.}
}
\description{
Before starting simulation, put data in canonical form, imputing
missing values and so forth as needed. Call compare_methods
repeatidly on the resulting treatment-permuted data.
}
