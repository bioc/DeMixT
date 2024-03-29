% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DeMixT_Preprocessing.R
\name{detect_suspicious_sample_by_hierarchical_clustering_2comp}
\alias{detect_suspicious_sample_by_hierarchical_clustering_2comp}
\alias{plot_sd}
\alias{plot_dim}
\title{detect_suspicious_sample_by_hierarchical_clustering_2comp}
\usage{
detect_suspicious_sample_by_hierarchical_clustering_2comp(
  count.matrix,
  normal.id,
  tumor.id
)

plot_sd(count.matrix, normal.id, tumor.id)

plot_dim(
  count.matrix,
  labels,
  legend.position = "bottomleft",
  legend.cex = 1.2
)
}
\arguments{
\item{count.matrix}{A matrix of raw expression count with \eqn{G} by \eqn{(My + M1)}, where \eqn{G} is the number
of genes, \eqn{My} is the number of mixed samples and \eqn{M1} is the number of normal samples. Row names are genes
 column names are sample ids.}

\item{normal.id}{A vector of normal sample ids}

\item{tumor.id}{A vector of tumor sample ids}

\item{legend.position}{Position of legend in the plot. Default is bottomleft.}

\item{legend.cex}{Character expansion factor relative to current par("cex"). Default = 1.2}
}
\value{
list object




}
\description{
Detect suspicious samples by a hierarchical clustering

This function is designed to evaluate the separation of tumor samples and normal samples in a PCA space. 
If some normal samples appear in the tumor-sample dominated cluster, these normal samples are likely to 
be tumor samples and they are supposed to be filtered out before downstream analysis. But for those tumor 
samples appearing in the normal-sample dominated cluster, we do not remove them since they may be the ones 
with low tumor purity.

Plot the standard deviation of log2 raw expression

Plot the distribution of tumor and normal samples in a 2D PCA space based on their expressions
}
