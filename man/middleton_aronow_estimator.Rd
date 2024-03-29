% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/design_based_estimators.R
\name{middleton_aronow_estimator}
\alias{middleton_aronow_estimator}
\title{Estimate ATE via Middleton & Aronow unbiased approach}
\usage{
middleton_aronow_estimator(
  formula,
  data = NULL,
  control_formula = NULL,
  aggregated = FALSE,
  include_block_estimates = FALSE
)
}
\arguments{
\item{formula}{Formula for outcome and treatment and nesting.  If
NULL, data is assumed to be in canonical form (see vignette for
further discussion).}
}
\value{
tibble of estimates using different varieties of the
  methods described in the paper.
}
\description{
This follows Middleton, J. A., & Aronow, P. M. (2015). Unbiased
Estimation of the Average Treatment Effect in Cluster-Randomized
Experiments. 6(1–2), 312. doi: 10.1515/spp-2013-0002
}
\details{
They introduce a Raj Difference estimator that is an extension of a
Horvitz-Thompson estimator where, in effect, we divide the sum of
the outcomes by a fixed constant rather than the realized sample
size to avoid biasing our estimate.

Note that the control_formula is ignored--these estimators do not
adjust for additioanl controls at this time.
}
