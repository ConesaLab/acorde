#' @import dplyr
#' @import purrr
#' @import tibble
#' @import SingleCellExperiment
#' @importFrom Rdpack reprompt



#' @title Reorder cell type-specific expression matrix during co-expression simulation
#'
#' @description This function is used internally by \code{acorde} to perform
#' the shuffling of simulated features for an individual cell type, as part of
#' the co-expression simulation process. The function is called recursively by
#' \code{\link[acorde:simulate_coexpression]{simulate_coexpression()}} to
#' perform the simulation on a full scRNA-seq matrix.
#'
#' @param sim_data A count matrix with features as rows and cells as columns.
#' Feature IDs must be included in an additional column named \code{feature}.
#'
#' @param feature_ids A two-column \code{tibble} containing \code{top} and \code{bottom}
#' columns, each including the feature IDs of features to be used as highly or
#' lowly expressed when shuffling by the indicated expression pattern.
#'
#' @param group_pattern A logical vector, containing \code{TRUE} to indicate that
#' high expression in that cell type is desired and \code{FALSE} if the opposite.
#' The vector must be ordered as the cell types in \code{sim_data}.
#'
#' @param cluster_no An integer indicating the total number of co-expression
#' clusters that are being generated in the simulation
#'
#' @return An expression matrix, with the same characteristics as \code{sim_data},
#' and a number of features defined as the total amount of top/bottom features
#' selected divided by the number of clusters for which co-expression patterns
#' where supplied.
#'
shuffle_group_matrix <- function(sim_data, feature_ids, group_pattern, cluster_no){

  # select top and bottom features in group
  top <- dplyr::select(feature_ids, top) %>% unlist
  bottom <- dplyr::select(feature_ids, bottom) %>% unlist

  # random partitioning of features
  # top
  top.shuffle <- sample(length(top))
  top <- top[top.shuffle]
  top.list <- split(top, cut(seq(1, length(top)), breaks = cluster_no, labels = FALSE))
  # bottom
  bottom.shuffle <- sample(length(bottom))
  bottom <- bottom[bottom.shuffle]
  bottom.list <- split(bottom, cut(seq(1, length(bottom)), breaks = cluster_no, labels = FALSE))

  # bind features following pattern
  features_bound <- vector(mode = "list", length = length(group_pattern))
  features_bound[group_pattern] <- top.list
  features_bound[!group_pattern] <- bottom.list
  features_bound <- unlist(features_bound)

  # build expression matrix for group
  sim_data.mod <- sim_data %>%
    dplyr::filter(feature %in% features_bound) %>%
    tibble::column_to_rownames("feature")
  sim_data.mod <- sim_data.mod[features_bound,] %>% tibble::rownames_to_column("feature")

  return(sim_data.mod)
}


#' @title Simulate co-expression patterns on an already-simulated scRNA-seq matrix
#'
#' @description Using a existing simulated scRNA-seq matrix, this function
#' creates co-expression relationships between the features, following the
#' cell type-specific patterns of high/low expression supplied by the user.
#' In particular, the output of the \code{SymSim}
#' simulator is expected as input (see \code{SymSim} package documentation for details).
#'
#' @param sim_data A \code{SingleCellExperiment} object containing an already-simulated
#' scRNA-seq matrix.
#'
#' @param feature_no An integer indicating the number of high expression ("top")
#' and low expression ("bottom") features to be selected for co-expression
#' simulation. Note that the output matrix will contain \code{feature_no*2}
#' features in total.
#'
#' @param patterns A \code{data.frame} or \code{tibble} containing cell types as
#' columns (ordered as in the \code{colData} slot in \code{sim_data}) and
#' co-expression clusters as rows. For each co-expression cluster,
#' a logical vector indicating the desired expression pattern must be
#' provided, in a row-wise manner. Insert \code{TRUE} if
#' high expression in that cell type is desired, \code{FALSE} if the opposite.
#'
#' @return A \code{list}, containing two objects:
#'
#' \enumerate{
#'
#' \item \code{sim_matrix}: a \code{tibble} containing the same number of cells
#' as in \code{sim_data} in the columns and \code{feature_no*2} in the rows.
#' Feature IDs are defined in the \code{feature} column
#'
#' \item \code{sim_clusters}: a \code{list} with as many elements as simulated
#' clusters, where each element contains all feature IDs that were simulated to
#' follow the same co-expression pattern (that is, the clusters).
#'
#' }
#'
#' @export
simulate_coexpression <- function(sim_data,
                                  feature_no,
                                  patterns){

  ## DATA PREPARATION: CELL TYPE SPECIFIC MATRICES AND FEATURES ##

  # extract counts from SCE
  normcounts <- SingleCellExperiment::counts(sim_data) %>% as.data.frame
  # get cell ids in each cell type
  group.list <- sim_data$Cell %>% split(sim_data$Group)
  # extract cell type (group) expr matrices
  normcounts.list <- purrr::map(group.list,
                                ~(normcounts[, as.character(.)] %>%
                                    tibble::rownames_to_column("feature")))

  # rank features by mean expression in each cell type
  normcounts.list <- map(normcounts.list,
                                ~dplyr::mutate(., mean = rowMeans(.[,-1])) %>%
                                  dplyr::arrange(dplyr::desc(mean))
                         %>% dplyr::select(-mean))

  # select top and bottom feature IDs for each cell type

    # top
    top_features.list <- purrr::map(normcounts.list,
                                ~dplyr::select(., feature) %>%
                                  dplyr::rename(top = "feature") %>%
                                  dplyr::slice(., 1:feature_no) %>%
                                  tibble::as_tibble())

    # modify bottom feature no. to create range correctly
    feature_no.c <- nrow(normcounts) - (feature_no - 1)
    # bottom
    bottom_features.list <- purrr::map(normcounts.list,
                                       ~dplyr::select(., feature) %>%
                                         dplyr::rename(bottom = "feature") %>%
                                         dplyr::slice(., feature_no.c:nrow(normcounts)) %>%
                                         tibble::as_tibble())

    # create a two-col tibble with top/bottom features per group
    features.list <- purrr::map2(top_features.list, bottom_features.list,
                                 dplyr::bind_cols)


  ## USE SUPPLIED PATTERS TO SHUFFLE THE CELL TYPE MATRICES ##

  # match column names for patterns
  colnames(patterns) <- names(normcounts.list)

  # shuffle matrix for each cell type following cluster patterns
  # note that internal function shuffle_group_matrix() is used to perform
  # each individual shuffling operation
  expr.list <- purrr::pmap(list(normcounts.list, features.list, patterns),
                    ~shuffle_group_matrix(sim_data = ..1,
                                          feature_ids = ..2,
                                          group_pattern = ..3,
                                          cluster_no = nrow(patterns)))

  # join cell type matrices into a single expression matrix
  expr.list <- purrr::map(expr.list, select, -feature)
  coexpr.df <- dplyr::bind_cols(expr.list) %>% tibble::as_tibble()
  coexpr.df <- coexpr.df %>%
    dplyr::mutate(feature = paste0("Feature", seq(1, nrow(coexpr.df)))) %>%
    dplyr::relocate(feature, .before = Cell1)

  # generate feature ID vectors for co-expression clusters
  clusters <- split(coexpr.df$feature,
                    cut(seq(1, nrow(coexpr.df)),
                        breaks = nrow(patterns), labels = FALSE))

  # build a list with results
  coexpr_sim <- list(sim_matrix = coexpr.df,
                     sim_clusters = clusters)
  return(coexpr_sim)

}
