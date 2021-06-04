#' @title Tasic data set (post-downsampling)
#'
#' @description Mouse neural single-cell RNA-Seq data set containing
#' isoform-level expression obtained using a bulk long-read transcriptome and
#' single-cell short-reads obtained from Tasic et al. (see reference below).
#'
#' @format A tibble containing single-cell isoform counts for 241 cells and
#' 13452 isoforms, represented as columns and rows, respectively.
#' Isoform identifiers are defined the \code{transcript_id} column. The data
#' has been length-normalized, quality filtered and subject to cell type-level
#' downsampling, as detailed in Arzalluz-Luque et al. (see reference below).
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
"tasic_down"


#' @title Tasic metadata (post-downsampling)
#'
#' @description Metadata table provided from publicly available experimental
#' details from Tasic et al. (see references).
#'
#' @format A tibble. Contains 10 columns, including cell IDs (taken from
#' sequencing run IDs in the \code{run} column), cell types and subtypes assgined
#' by the authors (\code{cell_type} and \code{subtype} columns).
#'
#'
#' @references
#' Public single-cell RNA-Seq data (Smart-Seq2) from mouse visual cortex was
#' obtained from Tasic et al.:
#'
#' \insertRef{Tasic2016a}{acorde}
"metadata_down"


#' @title Gene-isoform associations for the mouse neural long read transcriptome
#'
#' @description Table associating long read-defined isoform IDs to their genes.
#'
#' @format A tibble containing the columns \code{transcript_id} and \code{gene_id}.
#' Isoform (i.e. transcript) identifiers correspond to those obtained during
#' long read-transcriptome definition (see Arzalluz Luque et al.)
#'
#' @references
#'
#' Long read transcriptome definition details can be found in the methods section of
#' Arzalluz-Luque et al.:
#'
#' \insertRef{Arzalluz-Luque2021}{acorde}
