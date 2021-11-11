#' @title Tasic et al. 2016 expression matrix
#'
#' @description Mouse neural single-cell RNA-Seq data set containing
#' isoform-level expression. Expression matrix obtained using a bulk long-read
#' transcriptome and single-cell short-reads obtained from Tasic et al. 2016
#' (see reference below).
#'
#' @format A tibble containing single-cell isoform counts for 1591 cells and
#' 16240 isoforms, represented as columns and rows, respectively.
#' Isoform identifiers are defined the \code{transcript_id} column. The data
#' has been length-normalized and quality filtered as detailed in
#' Arzalluz-Luque et al. (see reference below).
#'
#' @references
#'
#' Long read transcriptome definition details can be found in the methods section of
#' Arzalluz-Luque et al.:
#'
#' \insertRef{Arzalluz-Luque2021}{acorde}
#'
#' Public single-cell RNA-Seq data (Smart-Seq2) from mouse visual cortex was
#' obtained from Tasic et al.:
#'
#' \insertRef{Tasic2016a}{acorde}
"tasic"


#' @title Tasic et al. 2016 metadata
#'
#' @description Metadata table provided from publicly available experimental
#' details from Tasic et al. 2016 (see references).
#'
#' @format A tibble. Contains 10 columns, including cell IDs (taken from
#' sequencing run IDs in the \code{run} column), cell types and subtypes assgined
#' by the authors (\code{cell_type} and \code{subtype} columns).
#'
#' @references
#' Public single-cell RNA-Seq data (Smart-Seq2) from mouse visual cortex was
#' obtained from Tasic et al.:
#'
#' \insertRef{Tasic2016a}{acorde}
"metadata"


#' @title Gene-isoform associations for the mouse neural long read transcriptome
#'
#' @description Table associating long read-defined isoform IDs to their genes.
#'
#' @format A tibble containing the columns \code{transcript} and \code{gene}.
#' Isoform (i.e. transcript) identifiers correspond to those obtained during
#' long read-transcriptome definition (see Arzalluz Luque et al.)
#'
#' @references
#'
#' Long read transcriptome definition details can be found in the methods section of
#' Arzalluz-Luque et al.:
#'
#' \insertRef{Arzalluz-Luque2021}{acorde}
"gene_tr_ID"


#' @title Isoform Differential Expression analysis results
#'
#' @description This object contains a list of differentially expressed isoforms
#' obtained after x50 downsamplping runs combined with Differential Expression
#' testing using edgeR and DESeq2. Given potential variability across
#' DE method runs and to ensure the reproducibility of the results in Arzalluz-Luque
#' et al. (see refs.), we provide the original results used in the manuscript.
#'
#' @format A tibble containing a single column, named \code{transcript}.
#' This transcript list is used to filter the \code{\link{tasic}} expression
#' matrix to retain only consistently DE transcript isoforms.
#'
#' @references
#'
#' \insertRef{Arzalluz-Luque2021}{acorde}
"consensus_DE_set"

