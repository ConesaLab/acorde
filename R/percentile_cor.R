#' @import dplyr
#' @import purrr
#' @import tibble
#' @import stringr


#' @title Summarize cell type-level expression using percentiles
#'
#' @description For each of the isoforms in the single-cell data set,
#' the function calculates percentile values using cell type-level expression values,
#' providing an estimate of the expression distribution shown by a particular
#' isoform in each of the cell types.
#'
#' @param data A data.frame or tibble object including isoforms as rows and
#' cells as columns. Isoform IDs can be included as row names (data.frame)
#' or as an additional column (tibble).
#' @param id_table A data frame including two columns named \code{cell} and
#' \code{cell_type}, in which correspondence between cell ID and cell type should be
#' provided. The number of rows should be equal to the total number of
#' cell columns in \code{data}, and the order of the \code{cell} column should
#' match column (i.e. cell) order in \code{data}.
#' @param percentile_no Integer indicating the number of percentiles to generate
#' for each of the cell types. Should always be higher than 4 (quantiles)
#' and lower than 100 (percentiles). Defaults to 10.
#' @param isoform_col When a tibble is provided in \code{data}, a character value
#' indicating the name of the column in which isoform IDs are specified.
#'
#' @return A \code{\link[tibble]{tibble}} containing one column of percentile-
#' summarized expression values per input transcript
#' in \code{data}.
#'
#' @export
percentile_expr <- function(data, id_table, percentile_no = 10, isoform_col = NULL){

  # handle rownames and data type
  if(is.null(isoform_col) == TRUE){
    data <- data %>% as.data.frame %>% rownames_to_column("transcript")
  }

  # split cell IDs by cell type labels
  cells_split <- split(id_table$cell, id_table$cell_type)

  # check that percentile_no is between 4 and 100
  if(percentile_no < 4 | percentile_no > 100){
    warning("percentile_no is not between 4 (quantiles) and 100 (percentiles).")
  }

  # generate step for probabilities in quantile function
  step <- 1/percentile_no
  # calculate percentiles by cell type for each transcript
  percentile_list <- map(cells_split, ~(select(data, all_of(.)) %>%
                                          apply(1, quantile, seq(0, 1, step)) %>% as.data.frame))
  percentiles <- bind_rows(percentile_list)
  colnames(percentiles) <- data[[1]]

  return(percentiles)
}


#' @title Compute percentile correlations between a set of transcripts
#'
#' @description This function summarizes expression for each cell type
#' using percentiles and then calculates Pearson correlations between all possible
#' transcript pairs in the data set. Internally, the function function calls
#' \code{\link{percentile_expr}} to generate percentile-summarized expression, and
#' then \code{\link[stats]{cor}} to compute Pearson correlation
#'
#' @param data A data.frame or tibble object including isoforms as rows and
#' cells as columns. Isoform IDs can be included as row names (data.frame)
#' or as an additional column (tibble).
#' @param id_table A data frame including two columns named \code{cell} and
#' \code{cell_type}, in which correspondence between cell ID and cell type should be
#' provided. The number of rows should be equal to the total number of
#' cell columns in \code{data}, and the order of the \code{cell} column should
#' match column (i.e. cell) order in \code{data}.
#' @param percentile_no Integer indicating the number of percentiles to generate
#' for each of the cell types. Should always be higher than 4 (quantiles)
#' and lower than 100 (percentiles). Defaults to 10.
#' @param isoform_col When a tibble is provided in \code{data}, a character value
#' indicating the name of the column in which isoform IDs are specified.
#'
#' @return A correlation matrix including one column-row pair per isoform
#' included in \code{data}.
#'
#' @export
percentile_cor <- function(data, id_table, percentile_no = 10, isoform_col = NULL){

  # get percentile expression
  percentiles <- percentile_expr(data, id_table, percentile_no, isoform_col)

  # calculate correlations
  cors <- stats::cor(percentiles)

  return(cors)
}
