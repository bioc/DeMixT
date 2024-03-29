% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DeMixT_Preprocessing.R
\name{DeMixT_preprocessing}
\alias{DeMixT_preprocessing}
\title{DeMixT_preprocessing}
\usage{
DeMixT_preprocessing(
  count.matrix,
  normal.id,
  tumor.id,
  cutoff_normal_range = c(0.1, 1),
  cutoff_tumor_range = c(0, 2.5),
  cutoff_step = 0.2
)
}
\arguments{
\item{count.matrix}{A matrix of raw expression count with \eqn{G} by \eqn{(My + M1)}, where \eqn{G} is the number
of genes, \eqn{My} is the number of mixed samples and \eqn{M1} is the number of normal samples. Row names are genes
 column names are sample ids.}

\item{normal.id}{A vector of normal sample ids}

\item{tumor.id}{A vector of tumor sample ids}

\item{cutoff_normal_range}{A vector of two numeric values, indicating the lower and upper bounds of standard deviation of 
log2 count matrix from the normal samples to subset. Default is c(0.2, 0.6)}

\item{cutoff_tumor_range}{A vector of two numeric values, indicating the lower and upper bounds to search standard deviation of 
log2 count matrix from the normal samples to subset. Default is c(0.2, 0.6)}

\item{cutoff_step}{A scatter value indicating the step size of changing cutoff_normal_range and cutoff_tumor_range to find a 
suitable subset of count matrix for downstream analysis}
}
\value{
processed count matrix
}
\description{
DeMixT preprocessing in one go
}
