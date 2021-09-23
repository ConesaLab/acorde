#' @import dplyr
#' @import purrr
#' @import tibble
#' @importFrom Rdpack reprompt
#'
#'
#' @title Filter isoforms by maximum cell type-level proportion of zeros
#'
#' @description For convenience, \code{acorde} includes a function to assist
#' lenient filtering of isoforms based on the proportion of zero expression values
#' across cell types. Based on this criteria, a minimum number of cells must have
#' non-zero expression in at least one cell type.
#'
#' @param data A data.frame or tibble object including isoforms as rows and cells as columns.
#' Isoform IDs can be included as row names (data.frame) or as an additional column (tibble).
#'
#' @param id_table A data frame including two columns named \code{cell} and \code{cell_type},
#' in which correspondence between cell ID and cell type should be
#' provided. The number of rows should be equal to the total number of
#' cell columns in \code{data}, and the order of the \code{cell} column should
#' match column (i.e. cell) order in \code{data}.
#'
#' @param ct_proportion A numeric indicating the minimum proportion of cells with
#' non-zero expression that will be allowed per cell type. Isoforms with a non-zero
#' value proportion above the threshold in at least one cell type will be flagged
#' to be preserved. Defaults to 0.2 (i.e. 20\%).
#'
#' @param isoform_col When a tibble is provided in \code{data}, a character object
#' indicating the name of the column in which isoform IDs are specified.
#' Otherwise, isoform identifiers will be assumed to be defined as rownames,
#' and this argument will not need to be provided.
#'
#' @return A logical vector including one entry per isoform in \code{data}. Isoforms
#' meeting the sparsity criteria will have a value of \code{TRUE}, and otherwise
#' be labeled as \code{FALSE}. This logical vector can then be used to filter
#' isoforms, i.e. the rows in \code{data}.
#'
#' @export
detect_sparse <- function(data, id_table, ct_proportion = 0.2, isoform_col = NULL){

  # handle rownames
  if(is.null(rownames) == FALSE){
    data <- data %>% as.data.frame %>% column_to_rownames(isoform_col)
  }

  # test number of zeros in each cell type
  split <- split(data %>% t %>% as.data.frame, id_table$cell_type)
  test_zero <- map(split, ~(. > 0))

  # compare to cell type proportion threshold
  lgl <- map(test_zero, ~(colSums(.) >= nrow(.)*ct_proportion))
  allct_lgl <- bind_rows(lgl)

  # no. of non-zero values should be higher than proportion in at least one cell type
  final_lgl <- rowSums(allct_lgl) >= 1

  return(final_lgl)
}


#' @title Filter isoforms by gene-relative expression across cell types
#'
#' @description For convenience, \code{acorde} includes a function to assist
#' lenient filtering of isoforms based on the mean gene-relative expression of
#' the isoform across cell types. Isoforms that contribute marginally to total gene
#' expression in all cell types will be flagged as minor. Based on this criteria,
#' isoforms must have higher mean gene-relative expression than specified in at
#' least one of the cell types.
#'
#' @param data A data.frame or tibble object including isoforms as rows and cells as columns.
#' Isoform IDs can be included as row names (data.frame) or as an additional column (tibble).
#'
#' @param id_table A data frame including two columns named \code{cell} and \code{cell_type},
#' in which correspondence between cell ID and cell type should be
#' provided. The number of rows should be equal to the total number of
#' cell columns in \code{data}, and the order of the \code{cell} column should
#' match column (i.e. cell) order in \code{data}.
#'
#' @param gene_tr_table A data.frame or tibble containing two columns
#' named \code{transcript_id} and \code{gene_id}, indicating gene-isoform
#' correspondence.
#'
#' @param gene_expr_proportion A numeric value indicating the mean gene-relative
#' expression threshold that will be used to detect minor isoforms across cell
#' types. Defaults to 0.2 (i.e. mean 20\% of the gene's total expression in at
#' least one cell type).
#'
#' @param isoform_col When a tibble is provided in \code{data}, a character object
#' indicating the name of the column in which isoform IDs are specified.
#' Otherwise, isoform identifiers will be assumed to be defined as rownames,
#' and this argument will not need to be provided.
#'
#' @return A \code{tibble} containing two columns, the first one including
#' transcript IDs, and the second containing logical values specifying whether
#' the isoform was flagged as minor (considering the provided gene-relative
#' expression threshold).
#'
#' @export
detect_minor_isoforms <- function(data, id_table,
                                  gene_tr_table,
                                  gene_expr_proportion = 0.2,
                                  isoform_col = NULL){

  # handle rownames
  if(is.null(isoform_col) == TRUE){
    data <- data %>% tibble::rownames_to_column("transcript")
  }else{
    data <- data %>% dplyr::rename(transcript = isoform_col)
  }

  # create long-formatted matrix
  data_long <- data %>%
    tidyr::pivot_longer(-transcript,
                        names_to = "cell", values_to = "expression") %>%
    dplyr::left_join(cell_types, by = "cell")

  # add gene ID
  data_long <- data_long %>%
    dplyr::left_join(gene_tr_table, by = "transcript")

  # compute gene expression
  gene_expr <- data_long %>%
    dplyr::group_by(cell, gene) %>%
    dplyr::summarize(gene_expr = sum(expression))

  # add gene expression to transcript table
  data_long <- data_long %>%
    dplyr::left_join(gene_expr, by = c("gene", "cell"))

  # calculate relative expression of transcripts in each cell
  data_long <- data_long %>%
    dplyr::mutate(gene_prop = dplyr::if_else(gene_expr > 0,
                                             true = expression/gene_expr,
                                             false = 0))

  # calculate mean relative expression for each cell type
  data_byct <- data_long %>%
    dplyr::group_by(transcript, cell_type) %>%
    dplyr::summarize(mean_gene_prop = mean(gene_prop))

  # calculate maximum mean relative expression across all cell types
  data_maxprop <- data_byct %>%
    dplyr::group_by(transcript) %>%
    dplyr::summarize(max_gene_prop = max(mean_gene_prop))

  # compare to provided threshold and create output tibble
  # TRUE if the maximum mean proportion is less than gene_expr_proportion
  minor_df <- data_maxprop %>%
    dplyr::mutate(minor_isoform = max_gene_prop < gene_expr_proportion) %>%
    select(transcript, minor_isoform)

  return(minor_df)
}



