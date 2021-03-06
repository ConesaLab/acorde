% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{tasic}
\alias{tasic}
\title{Tasic et al. 2016 expression matrix}
\format{
A tibble containing single-cell isoform counts for 1591 cells and
16240 isoforms, represented as columns and rows, respectively.
Isoform identifiers are defined the \code{transcript_id} column. The data
has been length-normalized and quality filtered as detailed in
Arzalluz-Luque et al. (see reference below).
}
\usage{
tasic
}
\description{
Mouse neural single-cell RNA-Seq data set containing
isoform-level expression. Expression matrix obtained using a bulk long-read
transcriptome and single-cell short-reads obtained from Tasic et al. 2016
(see reference below).
}
\references{
Long read transcriptome definition details can be found in the methods section of
Arzalluz-Luque et al.:

\insertRef{Arzalluz-Luque2021}{acorde}

Public single-cell RNA-Seq data (Smart-Seq2) from mouse visual cortex was
obtained from Tasic et al.:

\insertRef{Tasic2016a}{acorde}
}
\keyword{datasets}
