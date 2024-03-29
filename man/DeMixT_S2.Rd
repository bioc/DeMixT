% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DeMixT_S2.R
\name{DeMixT_S2}
\alias{DeMixT_S2}
\title{Deconvolves expressions of each individual sample for unknown component}
\usage{
DeMixT_S2(
  data.Y,
  data.N1,
  data.N2 = NULL,
  givenpi,
  nbin = 50,
  nthread = parallel::detectCores() - 1
)
}
\arguments{
\item{data.Y}{A SummarizedExperiment object of expression data from mixed
tumor samples. It is a \eqn{G} by \eqn{My} matrix where \eqn{G} is the
number of genes and \eqn{My} is the number of mixed samples. Samples with
the same tissue type should be placed together in columns.}

\item{data.N1}{A SummarizedExperiment object of expression data from
reference component 1 (e.g., normal). It is a \eqn{G} by \eqn{M1} matrix
where \eqn{G} is the number of genes and \eqn{M1} is the number of samples
for component 1.}

\item{data.N2}{A SummarizedExperiment object of expression data from
additional reference samples. It is a \eqn{G} by \eqn{M2} matrix where
\eqn{G} is the number of genes and \eqn{M2} is the number of samples for
component 2. Component 2 is needed only for running a three-component model.}

\item{givenpi}{A vector of proportions for all mixed tumor samples. In
two-component analysis, it gives the proportions of the unknown reference
component, and in three-component analysis, it gives the proportions for the
two known components.}

\item{nbin}{Number of bins used in numerical integration for computing
complete likelihood. A larger value increases accuracy in estimation but
increases the running time, especially in a three-component deconvolution
problem. The default is 50.}

\item{nthread}{The number of threads used for deconvolution when OpenMP is
available in the system. The default is the number of whole threads minus
one. In our no-OpenMP version, it is set to 1.}
}
\value{
\item{decovExprT}{A matrix of deconvolved expression profiles
corresponding to T-component in mixed samples for a given subset of genes.
Each row corresponds to one gene and each column corresponds to one sample.}
\item{decovExprN1}{A matrix of deconvolved expression profiles corresponding
to N1-component in mixed samples for a given subset of genes. Each row
corresponds to one gene and each column corresponds to one sample.}
\item{decovExprN2}{A matrix of deconvolved expression profiles corresponding
to N2-component in mixed samples for a given subset of genes in a
three-component setting. Each row corresponds to one gene and each column
corresponds to one sample.} \item{decovMu}{A matrix of estimated \eqn{Mu} of
log2-normal distribution for both known (\eqn{MuN1, MuN2}) and unknown
component (\eqn{MuT}). Each row corresponds to one gene.}
\item{decovSigma}{Estimated \eqn{Sigma} of log2-normal distribution for both
known (\eqn{SigmaN1, SigmaN2}) and unknown component (\eqn{SigmaT}). Each
row corresponds to one gene.}
}
\description{
This function is designed to estimate the deconvolved expressions of
individual mixed tumor samples for unknown component for each gene.
}
\examples{

# Example 1: two-component deconvolution given proportions 
  data(test.data.2comp)
  givenpi <- c(t(as.matrix(test.data.2comp$pi[-2,])))
  res.S2 <- DeMixT_S2(data.Y = test.data.2comp$data.Y, 
                      data.N1 = test.data.2comp$data.N1,
                      data.N2 = NULL, 
                      givenpi = givenpi, 
                      nbin = 50)
#                  
# Example 2: three-component deconvolution given proportions 
# data(test.data.3comp)
# givenpi = c(t(test.data.3comp$pi[-3,])) 
# res <- DeMixT_S2(data.Y = test.data.3comp$data.Y, 
#                  data.N1 = test.data.3comp$data.N1,
#                  data.N2 = test.data.3comp$data.N2, 
#                  givenpi = givenpi, 
#                  nbin = 50)


}
\references{
Wang Z, Cao S, Morris J S, et al. Transcriptome Deconvolution of
Heterogeneous Tumor Samples with Immune Infiltration. iScience, 2018, 9:
451-460.
}
\seealso{
http://bioinformatics.mdanderson.org/main/DeMixT
}
\author{
Zeya Wang, Wenyi Wang
}
\keyword{DeMixT_S2}
