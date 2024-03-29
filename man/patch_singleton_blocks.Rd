% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patch_data_set.R
\name{patch_singleton_blocks}
\alias{patch_singleton_blocks}
\title{Drop or combine all-tx or all-co blocks}
\usage{
patch_singleton_blocks(
  formula = NULL,
  data,
  drop_data = TRUE,
  pool_clusters = TRUE,
  warn_missing = TRUE
)
}
\arguments{
\item{drop_data}{Drop the troublesome blocks if TRUE, pool them if
FALSE.}

\item{pool_clusters}{If pooling blocks rather than dropping them,
also pool clusters in 100% tx or 100% co blocks into single
cluster.}

\item{warn_missing}{Say something if anything happens.}
}
\value{
data with the blockID column modified or troublesome rows
  of data missing.
}
\description{
Given dataset with some blocks that are all treated or all control,
drop those blocks or replace the block ID with canonical new block
ID shared by all such blocks.
}
\details{
Missing values in the block ID are all considered singletons.

Also, depending on pool_clusters, pool the clusters in each of
these identified blocks into single clusters.
}
