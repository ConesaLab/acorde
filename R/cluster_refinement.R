#' @import dplyr
#' @import purrr
#' @import tibble



#' @title Remove isoforms with low correlation from clusters
#'
#' @description Refine the expression signal of a set of isoform clusters
#' by removing elements with low correlations with the rest of isoforms in the
#' cluster.
#'
#' @param cluster_list A list of character vectors, each containing the
#' identifiers of the isoforms in a cluster.
#'
#' @param contains_unclustered A logical, indicating whether the first element
#' in \code{cluster_list} corresponds to unassigned isoforms
#' (i.e. unclustered, see \code{\link{cluster_isoforms}}).
#'
#' @param cor_matrix A matrix including co-expression values (e.g.
#' correlations) between isoforms, where column and row names indicate isoform
#' IDs.
#'
#' @param min_cor A numeric value indicating the minimum correlation
#' (or co-expression value) allowed between isoforms of the same cluster.
#'
#' @param lowcor_threshold An integer indicating the number of correlation values
#' below \code{min_cor} allowed per isoform in the cluster. As a result, any
#' isoform with \code{correlation < min_cor} with \code{>= lowcor_threshold}
#' counterparts will be moved to the unclustered group.
#'
#' @param size_filter A logical indicating whether size filter of the correlation-filtered
#' filtered clusters should be performed.
#' Size filter removes clusters that contain too few elements after filtering
#' and moves their isoforms to the unclustered group.
#'
#' @param size_threshold An integer value indicating the minimum number of
#' isoforms required to preserve the cluster during size filtering.
#'
#' @return A list of clusters, where the first element will correspond to
#' unclustered isoforms and the rest to the correlation-filtered clusters.
#'
#' @export
filter_clusters <- function(cluster_list, contains_unclustered = TRUE,
                            cor_matrix, min_cor = 0.8, lowcor_threshold = 3,
                             size_filter = TRUE, size_threshold = 3){

  # if unclustered transcripts are contained in cluster_list, handle
  if(contains_unclustered == TRUE){
    unclustered <- cluster_list[[1]]
    cluster_list <- cluster_list[-1]
  }

  # single cluster filtering function multiple times across cluster_list
  filtered_list <- map(cluster_list, single_cluster_filter,
                       cor_matrix = cor_matrix, min_cor = min_cor,
                       lowcor_threshold = lowcor_threshold)

  # move discarded transcripts to unclustered group
  if(contains_unclustered == TRUE){
    # join to current unclustered group
    unclustered <- c(unclustered,
                     map2(filtered_list, cluster_list, ~(.y[!(.y %in% .x)])) %>% unlist)

  }else{
    # if unclustered not provided, simply generate new unclustered group
    unclustered <- map2(filtered_list, cluster_list, ~(.y[!(.y %in% .x)]))
  }

  # filter out small clusters (by threshold) and move to unclustered
  if(size_filter == TRUE){
    unclustered <- c(unclustered,
                     filtered_list[map_int(filtered_list, length) <= size_threshold] %>%
                       unlist) %>% unlist
    unclustered <- list("0" = unclustered)

    filtered_list <- filtered_list[map_int(filtered_list, length) > size_threshold]
    names(filtered_list) <- seq(1, length(filtered_list))
  }

  return(c(unclustered, filtered_list))

}



#' @describeIn filter_clusters applies the filtering parameters defined in
#' \code{filter_clusters} to one cluster. This function is called internally and
#' used to iterate and filter all clusters supplied via the \code{cluster_list}
#' argument.
single_cluster_filter <- function(cluster,
                                  cor_matrix,
                                  min_cor,
                                  lowcor_threshold){

  # select correlations for transcripts in cluster
  clust_cors <- cor_matrix[cluster, cluster] %>% as_tibble
  # count no. of low correlations
  tr_filt <- map_lgl(clust_cors, ~(sum(. < min_cor) < lowcor_threshold))
  # filter the cluster
  filtered_cluster <- cluster[tr_filt]
  return(filtered_cluster)
}



##### FUNCTION TO JOIN UNCLUSTERED TRANSCRIPTS TO A SET OF EXTANT CLUSTERS: EXPAND #####

#' @export
expand_clusters <- function(data, isoform_col = NULL, id_table,
                            cluster_list, unclustered,
                            force_expand = TRUE, expand_threshold = NULL,
                            method = c("percentile", "pearson", "spearman",
                                       "rho", "zi_kendall"),
                            percentile_no = 10){

  message(paste("Expanding clusters using ", method, "correlation..."))

  if(method == "percentile"){

    # get percentile expression
    percentiles <- percentile_expr(data, id_table, percentile_no = percentile_no,
                                   isoform_col = isoform_col)

    # metatranscripts of clusters: compute mean-summarized percentile expression per transcript
    metatranscripts <- map(cluster_list,
                           ~(percentiles[,.] %>% as_tibble %>% rowMeans %>%
                               enframe(value = "mean_percentile")))
    metatr.df <- map(metatranscripts, select, mean_percentile) %>% bind_cols

    # calculate correlation between unclustered and metatranscripts
    unclust_percentiles <- percentiles[,unclustered]
    unclust_cor <- stats::cor(unclust_percentiles, metatr.df)

  }else if(method != "percentile"){

    if(is.null(isoform_col = TRUE)){
      data <- data %>% as.data.frame %>% rownames_to_column("transcript_id")
    }

    # get metatranscripts (mean expression of clusters)
    metatr.df <- map(cluster_list,
                     ~(filter(data, transcript_id %in% .) %>% column_to_rownames("transcript_id") %>%
                         colMeans %>% enframe(value = "cluster_mean", name = NULL))) %>% bind_cols

    # compute correlation of unclustered with metatranscripts
    unclust_expr <- data %>% filter(transcript_id %in% unclustered) %>% column_to_rownames("transcript_id") %>%
      as.matrix %>% t() %>% as.data.frame

    if(method == "rho"){
      mat <- bind_cols(metatr.df, unclust_expr) %>% as.matrix()

      # use dismay function to compute rho
      unclust_cor <- dismay::dismay(mat, metric = "rho_p", select = colnames(mat))
      unclust_cor <- unclust_cor[colnames(unclust_expr), colnames(metatr.df)]

    }else if(method == "zi_kendall"){
      mat <- bind_cols(metatr.df, unclust_expr) %>% as.matrix()

      # use dismay function to compute rho
      unclust_cor <- dismay::dismay(mat, metric = "zi_kendall")
      unclust_cor <- unclust_cor[colnames(unclust_expr), colnames(metatr.df)]

    }else{
      unclust_cor <- stats::cor(unclust_expr, metatr.df, method = method)
    }
  }


  # evaluate force_expand and leave unclustered if correlations are too low
  if(force_expand == FALSE){
    highcor <- apply(unclust_cor, 1, function(x) sum(x >= expand_threshold))
    if(sum(highcor) == 0) stop("No correlations above threshold found for unclustered features.")
    # filter correlation matrix to remove isoforms with no high correlation
    unclust_cor <- unclust_cor[names(highcor[highcor >= 1]),]
    # new unclustered group
    new_unclust <- names(highcor[highcor == 0])
  }

  # find maximum correlation and assign: transcripts go to cluster with maximum metatranscript correlation
  # if force_expand = FALSE, correlation matrix has been previously filtered
  assign_unclust <- apply(unclust_cor, 1, which.max) %>% split(names(.), .)

  # skip clusters with no assigned unclustered transcripts
  skip <- which(!(seq(1,length(cluster_list)) %in% (names(assign_unclust) %>% as.integer)))

  # assign unclustered
  if(length(skip) == 0){
    cluster_list.expanded <- map2(cluster_list, assign_unclust,
                                  ~(c(.x, .y)))
    cluster_list.nonexpanded <- NULL
  } else {
    cluster_list.expanded <- map2(cluster_list[-c(skip)], assign_unclust,
                                  ~(c(.x, .y)))
    # recover skipped clusters, i.e. clusters that have not been expanded
    cluster_list.nonexpanded <- cluster_list[skip]
  }

  # join and reorder by name (i.e. by original cluster number)
  cluster_list.final <- c(cluster_list.expanded, cluster_list.nonexpanded)

  order <- names(cluster_list.final) %>% as.integer %>% sort() %>% as.character
  cluster_list.final.ord <- cluster_list.final[order]

  # return depending on force_expand
  if(force_expand == FALSE){
    return(list(unclustered = new_unclust, expanded = cluster_list.final.ord))
  } else {
    return(cluster_list.final.ord)
  }
}



#### FUNCTION TO FILTER CLUSTERS BY DS AND SPLICING COORDINATION #####

#' @export
filter_coDIU <- function(cluster_list, gene_tr_table){

  message(paste("Total no. of clusters:", length(cluster_list), sep = " "))
  message(paste("Total isoforms in clusters:", unlist(cluster_list) %>% length), sep = " ")

  # filter out transcripts from genes with only one isoform left

  # list of all clustered transcripts
  clustered_tr <- unlist(cluster_list) %>% unname
  # split/group transcripts by gene IDs
  clustered_g <- split(clustered_tr, gene_tr_table[match(clustered_tr, gene_tr_table$transcript_id),]$gene_id)
  # count no. of transcripts per gene in the clusters
  ntr <- map_int(clustered_g, length)
  # find genes with > 1 isoform clustered
  gmulti <- ntr[ntr > 1]
  # find transcripts from genes with > 1 isoform clustered
  multi_tr <- unlist(clustered_g[names(gmulti)])
  # filter clusters to keep only transcripts with multiple same-gene transcripts
  clusters_multi <- map(cluster_list, ~(.[. %in% multi_tr]))

  message(paste("Isoforms clustered after coordination filter:", unlist(clusters_multi) %>% length, sep = " "))

  # filter out transcripts from genes with all isoforms in same cluster

  # convert clusters to gene IDs
  clusters_multi.gene <- map(clusters_multi, ~(gene_tr_table[match(., gene_tr_table$transcript_id),]$gene_id))

  # find no. of clusters where each gene has isoforms
  gene_distribution <- map(clusters_multi.gene, ~(names(gmulti) %in% .)) %>% bind_rows %>% t
  colnames(gene_distribution) <- names(gmulti)
  gene_distribution <- colSums(gene_distribution)
  # find differentially spliced genes (i.e. isoforms in more than one cluster)
  genes_ds <- gene_distribution[gene_distribution > 1]

  # filter clusters to only keep transcripts from ds genes
  tr_ds.idx <- map(clusters_multi.gene, ~(which(. %in% names(genes_ds))))
  clusters_ds <- map2(clusters_multi, tr_ds.idx, ~(.x[.y]))

  message(paste("Isoforms clustered after differential splicing filter:", unlist(clusters_ds) %>% length, sep = " "))

  return(clusters_ds)
}
