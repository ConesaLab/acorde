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
