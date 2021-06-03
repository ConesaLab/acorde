#' @import dplyr
#' @import purrr
#' @import tibble


#' @title Initial clustering of isoforms using dynamicTreeCut
#'
#' @description Create groups of isoforms with similar expression across cell
#' types. Based on co-expression values and dynamic
#' hierarchical clustering with the \code{dynamicTreeCut} package.
#' Ultimately, this function provides a custom wrapper for the
#' \code{\link[dynamicTreeCut]{cutreeHybrid}} function.
#'
#' @param cor_matrix A matrix including co-expression values (e.g.
#' correlations) between isoforms, where column and row names indicate isoform
#' IDs.
#'
#' @param deepSplit Either logical or integer int he range 0 to 4 (defaults to 3).
#' Provides control over sensitivity to cluster splitting (higher value creates
#' more clusters with a smaller number of elements).
#'
#' @param pamStage A logical indicating whether a second clustering step similar
#' to Partition Around Medioids is to be performed. If \code{pamStage = TRUE} is
#' set, no isoforms will remain unassigned to clusters. By default, we allow
#' unclustered elements by setting \code{pamStage = FALSE}.
#'
#' @param minClusterSize An integer defining minimum cluster size allowed
#' during dynamic clustering. Defaults to 20.
#'
#' @param ... Additional parameters supplied to
#' \code{\link[dynamicTreeCut]{cutreeHybrid}}.
#'
#' @details First, a previously obtained \strong{co-expression matrix} (for instance,
#' percentile correlations generated using \code{\link{percentile_cor}})
#' is transformed into a distance matrix for hierarchical clustering.
#' As the \code{acorde} pipeline intends to detect positively correlated isoforms,
#' negative correlation values are automatically discarded by replacing them
#' with zero values. Distance is computed as \emph{1 - |co-expression|}.
#'
#' Next, \strong{hierarchical clustering} is performed using the
#' \code{\link[stats]{hclust}} function. Using the output dendrogram, isoforms
#' are clustered via the \code{\link[dynamicTreeCut]{cutreeHybrid}} function
#' in the \code{dynamicTreeCut} package.
#'
#' @seealso For more details on parameters, refer to
#' \code{\link[dynamicTreeCut]{cutreeHybrid}}.
#'
#' @return A list of generated clusters, that is, of character vectors including
#' the identifiers of isoforms assigned to each of the clusters.
#'
#' @references
#' \insertRef{Langfelder2008}{acorde}
#'
#' \insertRef{Venables2002}{acorde}
#'
#' @export
cluster_isoforms <- function(cor_matrix, deepSplit = 3, pamStage = FALSE,
                             minClusterSize = 20,...){

  # hierarchical clustering on correlation-based distance
  message("Inferring dendrogram via hclust()...")
  dist <- stats::as.dist(1 - abs(cor_matrix))
  h <- stats::hclust(dist, method = "average")

  # create clusters dynamically from hclust dendrogram
  message("Creating clusters dynamically via cutreeHybrid()...")
  hdynamic <- dynamicTreeCut::cutreeHybrid(h, as.matrix(dist),
                           deepSplit = deepSplit, pamStage = pamStage,
                           minClusterSize = minClusterSize, ...)
  hclusters <- split(h$labels, hdynamic$labels)

  return(hclusters)
}



#' @title Merge redundant clusters by expression profile similarity
#'
#' @description Join clusters representing the same expression
#' pattern across cell types (redundant clusters). This function uses a
#' metaclustering system (see details) and user-defined similarity thresholds
#' that allows to control for the stringency of the merge process.
#'
#' @param data A data.frame or tibble object including isoforms as rows and
#' cells as columns. Isoform IDs can be included as row names (data.frame)
#' or as an additional column (tibble).
#'
#' @param isoform_col When a tibble is provided in \code{data}, a character value
#' indicating the name of the column in which isoform IDs are specified.
#'
#' @param id_table A data frame including two columns named \code{cell} and
#' \code{cell_type}, in which correspondence between cell ID and cell type should be
#' provided. The number of rows should be equal to the total number of
#' cell columns in \code{data}, and the order of the \code{cell} column should
#' match column (i.e. cell) order in \code{data}.
#'
#' @param cluster_list A list of character vectors, each containing the
#' identifiers of the isoforms in a cluster.
#'
#' @param percentile_no Integer indicating the number of percentiles that will
#' be used to summarized cell type expression via \code{\link{percentile_expr}}.
#' Should always be higher than 4 (quantiles) and lower than 100 (percentiles).
#' Defaults to 10.
#'
#' @param method Character indicating a co-expression method to use for merging
#' similar clusters. Should be one of \code{percentile, pearson, spearman,
#' zi_kendall, rho} (see details).
#'
#' @param dynamic A logical. If \code{TRUE}, merge will be performed
#' via dynamic hierarchical clustering. Defaults to \code{FALSE}.
#'
#' @param height_cutoff When \code{dynamic = FALSE}, a numeric value between
#' 0 and 1 to be supplied to \code{\link[stats]{cutree}} via the \code{h} argument.
#' Indicates the height where the created dendrogram tree should be cut to
#' generate groups of merged clusters.
#'
#' @param cutree_no An integer indicating the desired number of groups to
#' merge clusters into. Supplied to \code{\link[stats]{cutree}} via the
#' \code{k} argument. Only required when \code{dynamic = FALSE} and
#' \code{height_cutoff = NULL}.
#'
#' @param ... Additional arguments passed to \code{\link[dynamicTreeCut]{cutreeHybrid}}
#' (only when \code{dynamic = TRUE}).
#'
#' @details During the isoform clustering process, it is generally useful to
#' prioritize the reduction of within-cluster variability. This, however, can lead
#' to obtaining a large number of small, redundant clusters. To mitigate this
#' effect, \code{acorde} includes a step where clusters with high profile similarity
#' can be merged using the correlation between their \emph{metatranscripts}.
#' A cluster's metatranscript is calculated as the mean of the
#' \link[acorde:percentile_expr]{percentile-summarized}
#' expression of all of the isoforms in that cluster. Then, co-expression values
#' between metatranscripts are calculated and used to
#' generate a distance matrix to group cluster profiles by similarity, a process
#' that can be referred to as \emph{metaclustering}.
#'
#' By default, the \strong{metaclustering} proccess is done using traditional
#' hierarchical clustering via \code{\link[stats]{hclust}},
#' which requires the definition of either a height cutoff (\code{height_cutoff}
#' parameter) or a number of clusters to obtain (\code{cutree_no}).
#'
#' Available \strong{co-expression metrics} (selected via the \code{method}) include:
#' \enumerate{
#'    \item \code{percentile}: percentile correlations computed using
#'    \code{\link{percentile_cor}}.
#'
#'    \item \code{pearson}: Pearson correlation computed using
#'    \code{\link[stats]{cor}}.
#'
#'    \item \code{spearman}: Spearman correlation computed using
#'    \code{\link[stats]{cor}}.
#'
#'    \item \code{zi_kendall}: zero-inflated Kendall correlation computed
#'    using the \code{\link[dismay]{dismay}} function.
#'
#'    \item \code{rho}: rho proportionality metric computed using the
#'    \code{\link[dismay]{dismay}} function.
#' }
#'
#' Alternatively, users may choose to perform metatranscript clustering
#' dynamically using the \code{dynamicTreeCut} package, therefore setting \code{dynamic = TRUE}.
#' In this case, additional parameters will need to be supplied to the
#' \code{\link[dynamicTreeCut]{cutreeHybrid}} function via the \code{...} argument.
#' Note that \code{minClusterSize = 1} is set internally to allow clusters to
#' remain unmerged if no redundancies with the profiles of other clusters are
#' found.
#'
#' @references
#' \insertRef{Langfelder2008}{acorde}
#'
#' \insertRef{Venables2002}{acorde}
#'
#' \insertRef{Skinnider2019}{acorde}
#'
#' @export
merge_clusters <- function(data, isoform_col = NULL, id_table,
                           cluster_list,
                           percentile_no = 10,
                           dynamic = FALSE,
                           method = c("percentile", "pearson", "spearman",
                                      "rho", "zi_kendall"),
                           height_cutoff = 0.2,
                           cutree_no = NULL, ...){

  if(method == "percentile"){

    # get percentile expression
    percentiles <- percentile_expr(data, ids_to_type, percentile_no = percentile_no,
                                   isoform_col = isoform_col)

    # metatranscripts of clusters: compute mean-summarized percentile expression per transcript
    metatranscripts <- map(cluster_list,
                           ~(percentiles[,.] %>% as_tibble %>% rowMeans %>%
                               enframe(value = "mean_percentile")))
    metatr.df <- map(metatranscripts, select, mean_percentile) %>% bind_cols
    colnames(metatr.df) <- seq(1, length(cluster_list))

    # get correlation between metatranscripts of all clusters
    cors.meta <- stats::cor(metatr.df)
    # discard negative correlations
    cors.meta[cors.meta < 0] <- 0

  } else if(method != "percentile"){

    # get metatranscripts (mean expression of clusters)
    metatr.df <- map(cluster_list,
                     ~(data[.,] %>% colMeans %>%
                         enframe(value = "cluster_mean", name = NULL))) %>%
      bind_cols
    colnames(metatr.df) <- seq(1, length(cluster_list))

    if(method == "rho"){
      # use dismay function to calculate correlation between metatranscripts
      cors.meta <- dismay::dismay(metatr.df %>% as.matrix, metric = "rho_p",
                                  select = colnames(metatr.df))

    }else if(method == "zi_kendall"){
      # use dismay function to calculate correlation between metatranscripts
      cors.meta <- dismay::dismay(metatr.df %>% as.matrix, metric = "zi_kendall")

    }else {
      # get correlation between metatranscripts of all clusters
      cors.meta <- stats::cor(metatr.df, method = method)
      # discard negative correlations
      cors.meta[cors.meta < 0] <- 0
    }
  }

  # clustering of the metatranscritps (i.e. the cluster profiles) by correlation
  dist.meta <- stats::as.dist(1 - cors.meta)
  h.meta <- stats::hclust(dist.meta, method = "complete")
  # create groups from the dendrogram
  if(dynamic == TRUE){
    hcut <- dynamicTreeCut::cutreeHybrid(h.meta, as.matrix(dist.meta),
                                         minClusterSize = 1, ...)
    hcut <- hcut$labels

  } else if(dynamic == FALSE){
    hcut <- stats::cutree(h.meta, h = height_cutoff, k = cutree_no)
  }

  # merge groups of similar clusters
  merge <- split(h.meta$labels, hcut) %>% map(as.integer)
  clusters_merged <- map(merge, ~(cluster_list[.] %>% unlist))

  # return a list of merged groups as well as the clusters after merging
  return(list(merged_groups = merge, clusters = clusters_merged))
}
