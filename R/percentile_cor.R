##### FUNCTIONS TO COMPUTE PERCENTILE CORRELATIONS BETWEEN ISOFORMS #####
# data: a matrix or data frame with transcript IDs in rownames (or, in the case of a data frame, in the first column
# i.e. see rownames_to_column() function from the tibble package)
# ids_to_type: a data frame with cell ids in the first column and cell type labels in the second column

#' @export
percentile_expr <- function(data, ids_to_type, percentile_no){

  require(tidyverse)

  # handle rownames and data type
  if(str_detect(colnames(data), "transcript") %>% sum == 0 || is.matrix(data) == TRUE){
    data <- data %>% as.data.frame %>% rownames_to_column("transcript")
  }

  # split cell IDs by cell type labels
  cells_split <- split(ids_to_type[[1]], ids_to_type[[2]])

  # check that percentile_no is between 4 and 100
  if(percentile_no < 4 | percentile_no > 100){
    warning("percentile_no is not between 4 (quantiles) and 100 (percentiles).")
  }

  # generate step for probabilities in quantile function
  step <- 1/percentile_no
  # calculate percentiles by cell type for each transcript
  percentile_list <- purrr::map(cells_split, ~(dplyr::select(data, all_of(.)) %>%
                                          apply(1, quantile, seq(0, 1, step)) %>% as.data.frame))
  percentiles <- bind_rows(percentile_list)
  colnames(percentiles) <- data[[1]]

  return(percentiles)
}


#' @export
percentile_cor <- function(data, ids_to_type, percentile_no = 10){

  # get percentile expression
  percentiles <- percentile_expr(data, ids_to_type, percentile_no)

  # calculate correlations
  cors <- cor(percentiles)

  return(cors)
}
