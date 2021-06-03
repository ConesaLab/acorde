#' @title Tasic data set
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
