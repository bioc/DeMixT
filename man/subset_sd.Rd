% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DeMixT_Preprocessing.R
\name{subset_sd}
\alias{subset_sd}
\title{subset_sd}
\usage{
subset_sd(
  count.matrix,
  normal.id,
  tumor.id,
  cutoff_normal = c(0.1, 0.6),
  cutoff_tumor = c(0.2, 0.8)
)
}
\arguments{
\item{count.matrix}{A matrix of raw expression count with \eqn{G} by \eqn{(My + M1)}, where \eqn{G} is the number
of genes, \eqn{My} is the number of mixed samples and \eqn{M1} is the number of normal samples. Row names are genes
 column names are sample ids.}

\item{normal.id}{A vector of normal sample ids}

\item{tumor.id}{A vector of tumor sample ids}

\item{cutoff_normal}{A vector of two numeric values, indicating the lower and upper bounds of standard deviation of 
log2 count matrix from the normal samples to subset. Default is c(0.1, 0.6)}

\item{cutoff_tumor}{A vector of two numeric values, indicating the lower and upper bounds of standard deviation of 
log2 count matrix from the tumor samples to subset. Default is c(0.2, 0.8)}
}
\value{
A subset of the count matrix
}
\description{
Subset a count matrix given the the ranges of the standard deviations of the 
log2 expressions from the tumor and normal samples
}
