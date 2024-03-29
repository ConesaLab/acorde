% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{symsim_sce}
\alias{symsim_sce}
\title{SymSim simulated dataset}
\format{
A \code{SingleCellExperiment} object containing two assays, \code{counts}
and \code{logcounts}, where log-transformed counts are included for PCA
visualization, among other functions. The \code{colData} slot contains a
two-column dataframe with \code{Cell} and \code{Group} identifiers. PCA
is already calculated and stored in the \code{reducedDimNames} slot.
}
\usage{
symsim_sce
}
\description{
A \code{SingleCellExperiment} object containing the results of
simulating data with SymSim, as described in the co-expression simulation
vignette in \code{acorde}. This dataset contains 1000 cells divided equally
into 8 cell types, and 8000 genes/features.
}
\keyword{datasets}
