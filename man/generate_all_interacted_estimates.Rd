% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_interacted_estimators.R
\name{generate_all_interacted_estimates}
\alias{generate_all_interacted_estimates}
\title{Average block-level estimates to get overall ATE}
\usage{
generate_all_interacted_estimates(
  fitModel,
  data,
  use_full_vcov = FALSE,
  method = "LR-FI-CRVE",
  aggregated = FALSE
)
}
\arguments{
\item{fitModel}{The model that has the interacted estimates in it.}

\item{use_full_vcov}{TRUE/FALSE. When calculating standard errors,
should the full variance-covariance matrix of the block-level
estimates be used, or just the diagonal of standard errors?}

\item{method}{Prefix of the method.  Will add the weighting to stem.}
}
\description{
Also calculate associated standard error.
}