##### FUNCTION TO CLUSTER ISOFORMS BASED ON CORRELATION #####

# cor_matrix: matrix of correlation values with column names and rownames indicating the transcript IDs
# ...: can be used to change parameters of cutreeHybrid() function. Three set by default can also be changed.
cluster_isoforms <- function(cor_matrix, deepSplit = 3, pamStage = FALSE, minClusterSize = 20,...){
  
  require(dynamicTreeCut)
  
  # hierarchical clustering on correlation-based distance
  message("Inferring dendrogram via hclust()...")
  dist <- as.dist(1 - abs(cor_matrix))
  h <- hclust(dist, method = "average")
  
  # create clusters dynamically from hclust dendrogram
  message("Creating clusters dynamically via cutreeHybrid()...")
  hdynamic <- cutreeHybrid(h, as.matrix(dist), 
                           deepSplit = deepSplit, pamStage = pamStage, minClusterSize = minClusterSize, ...)
  hclusters <- split(h$labels, hdynamic$labels)
  
  return(hclusters)
}

#### FUNCTION TO FILTER ONE CLUSTER AND GET A CLEANER SIGNAL ####
single_cluster_filter <- function(cluster, cor_matrix, min_cor, lowcor_threshold){
  
  require(tidyverse)
  
  # select correlations for transcripts in cluster
  clust_cors <- cor_matrix[cluster, cluster] %>% as_tibble
  # count no. of low correlations
  tr_filt <- purrr::map_lgl(clust_cors, ~(sum(. < min_cor) < lowcor_threshold))
  # filter the cluster
  filtered_cluster <- cluster[tr_filt]
  return(filtered_cluster)
}

#### FUNCTION TO QC FILTER ALL CLUSTERS AT ONCE BY THE SPECIFIED PARAMETERS ####
# main filtering function that calls single_cluster_filter()
# specified correlation value (min_cor) must be satisifed for each isoform vs all in the cluster 
# allowing lowcor_threshold exceptions in the rule
# length filter eliminates small clusters and moves members to unclustered list supplied/created
filter_clusterQC <- function(cluster_list, cor_matrix, min_cor = 0.8, lowcor_threshold = 3,
                            contains_unclustered = TRUE, length_filter = TRUE, length_threshold = 3){
  
  require(tidyverse)
  
  # if unclustered transcripts are contained in cluster_list, handle
  if(contains_unclustered == TRUE){
    unclustered <- cluster_list[[1]]
    cluster_list <- cluster_list[-1]
  }
  
  # single cluster filtering function multiple times across cluster_list
  filtered_list <- purrr::map(cluster_list, single_cluster_filter, 
                       cor_matrix = cor_matrix, min_cor = min_cor, lowcor_threshold = lowcor_threshold)
  
  # move discarded transcripts to unclustered group
  if(contains_unclustered == TRUE){
    # join to current unclustered group
    unclustered <- c(unclustered,
                     purrr::map2(filtered_list, cluster_list, ~(.y[!(.y %in% .x)])) %>% unlist)
  
  }else{
    # if unclustered not provided, simply generate new unclustered group
    unclustered <- purrr::map2(filtered_list, cluster_list, ~(.y[!(.y %in% .x)]))
  }

  # filter out small clusters (by threshold) and move to unclustered
  if(length_filter == TRUE){
    unclustered <- c(unclustered, 
                     filtered_list[purrr::map_int(filtered_list, length) <= length_threshold] %>% unlist) %>% unlist
    unclustered <- list("0" = unclustered)
    
    filtered_list <- filtered_list[purrr::map_int(filtered_list, length) > length_threshold]
    names(filtered_list) <- seq(1, length(filtered_list))
  }
  
  return(c(unclustered, filtered_list))
  
}


##### FUNCTION TO JOIN UNCLUSTERED TRANSCRIPTS TO A SET OF EXTANT CLUSTERS: EXPAND #####
expand_clusters <- function(data, cluster_list, unclustered, ids_to_type, 
                            force_expand = TRUE, expand_threshold = NULL, 
                            method = c("percentile", "pearson", "spearman", "rho", "zi_kendall"), 
                            percentile_no = 10){
  
  require(tidyverse)
  require(dismay)
  
  message(paste("Expanding clusters using ", method, "correlation..."))
  
  if(method == "percentile"){
    
    # get percentile expression
    percentiles <- percentile_expr(data, ids_to_type, percentile_no = percentile_no)
    
    # metatranscripts of clusters: compute mean-summarized percentile expression per transcript
    metatranscripts <- purrr::map(cluster_list, 
                           ~(percentiles[,.] %>% as_tibble %>% rowMeans %>% enframe(value = "mean_percentile")))
    metatr.df <- purrr::map(metatranscripts, dplyr::select, mean_percentile) %>% bind_cols
    
    # calculate correlation between unclustered and metatranscripts
    unclust_percentiles <- percentiles[,unclustered]
    unclust_cor <- cor(unclust_percentiles, metatr.df)
    
  }else if(method != "percentile"){
    
    if(str_detect(colnames(data), "transcript") %>% sum == 0){
      data <- data %>% as.data.frame %>% rownames_to_column("transcript_id")
    }
    
    # get metatranscripts (mean expression of clusters)
    metatr.df <- purrr::map(cluster_list, 
                      ~(filter(data, transcript_id %in% .) %>% column_to_rownames("transcript_id") %>% 
                        colMeans %>% enframe(value = "cluster_mean", name = NULL))) %>% bind_cols
    
    # compute correlation of unclustered with metatranscripts
    unclust_expr <- data %>% filter(transcript_id %in% unclustered) %>% column_to_rownames("transcript_id") %>% 
      as.matrix %>% t() %>% as.data.frame
    
    if(method == "rho"){
      mat <- bind_cols(metatr.df, unclust_expr) %>% as.matrix()
      
      # use dismay function to compute rho    
      unclust_cor <- dismay(mat, metric = "rho_p", select = colnames(mat))
      unclust_cor <- unclust_cor[colnames(unclust_expr), colnames(metatr.df)]
      
    }else if(method == "zi_kendall"){
      mat <- bind_cols(metatr.df, unclust_expr) %>% as.matrix()
      
      # use dismay function to compute rho    
      unclust_cor <- dismay(mat, metric = "zi_kendall")
      unclust_cor <- unclust_cor[colnames(unclust_expr), colnames(metatr.df)]
      
    }else{
      unclust_cor <- cor(unclust_expr, metatr.df, method = method)
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
    cluster_list.expanded <- purrr::map2(cluster_list, assign_unclust, 
                                  ~(c(.x, .y)))
    cluster_list.nonexpanded <- NULL
  } else {
    cluster_list.expanded <- purrr::map2(cluster_list[-c(skip)], assign_unclust, 
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


##### FUNCTION TO MERGE SIMILAR CLUSTERS ####
# heigth_cutoff: value used in the cutree() function to make cluster merge. 
# *Try several cutoffs on data, default is not optimal for all datasets.
# dynamic merge usingn cutTreeHybrid() also implemented
merge_clusters <- function(data, cluster_list, ids_to_type, height_cutoff = 0.2, cutree_no = NULL, 
                           percentile_no = 10, dynamic = FALSE, 
                           method = c("percentile", "pearson", "spearman", "rho", "zi_kendall"), ...){
  
  require(tidyverse)
  
  if(method == "percentile"){
    
    # get percentile expression
    percentiles <- percentile_expr(data, ids_to_type, percentile_no = percentile_no)
    
    # metatranscripts of clusters: compute mean-summarized percentile expression per transcript
    metatranscripts <- purrr::map(cluster_list, 
                           ~(percentiles[,.] %>% as_tibble %>% rowMeans %>% enframe(value = "mean_percentile")))
    metatr.df <- purrr::map(metatranscripts, dplyr::select, mean_percentile) %>% bind_cols
    colnames(metatr.df) <- seq(1, length(cluster_list))
    
    # get correlation between metatranscripts of all clusters
    cors.meta <- cor(metatr.df)
    # discard negative correlations
    cors.meta[cors.meta < 0] <- 0
    
  } else if(method != "percentile"){
    
    # get metatranscripts (mean expression of clusters)
    metatr.df <- purrr::map(cluster_list, 
                     ~(data[.,] %>% colMeans %>% enframe(value = "cluster_mean", name = NULL))) %>% bind_cols
    colnames(metatr.df) <- seq(1, length(cluster_list))
    
    if(method == "rho"){
      # use dismay function to calculate correlation between metatranscripts
      cors.meta <- dismay(metatr.df %>% as.matrix, metric = "rho_p", select = colnames(metatr.df))
      
    }else if(method == "zi_kendall"){
      # use dismay function to calculate correlation between metatranscripts
      cors.meta <- dismay(metatr.df %>% as.matrix, metric = "zi_kendall")
      
    }else {
      # get correlation between metatranscripts of all clusters
      cors.meta <- cor(metatr.df, method = method)
      # discard negative correlations
      cors.meta[cors.meta < 0] <- 0
    }
  }
  
  # clustering of the metatranscritps (i.e. the cluster profiles) by correlation
  dist.meta <- as.dist(1 - cors.meta)
  h.meta <- hclust(dist.meta, method = "complete")
  # create groups from the dendrogram
  if(dynamic == TRUE){
    require(dynamicTreeCut)
    hcut <- cutreeHybrid(h.meta, as.matrix(dist.meta), minClusterSize = 1, ...)
    hcut <- hcut$labels
  } else if(dynamic == FALSE){
    hcut <- cutree(h.meta, h = height_cutoff, k = cutree_no)
  }
  
  # merge groups of similar clusters
  merge <- split(h.meta$labels, hcut) %>% purrr::map(as.integer)
  clusters_merged <- purrr::map(merge, ~(cluster_list[.] %>% unlist))
  
  # return a list of merged groups as well as the clusters after merging
  return(list(merged_groups = merge, clusters = clusters_merged))
}



#### FUNCTION TO FILTER CLUSTERS BY DS AND SPLICING COORDINATION #####
filter_coDS <- function(cluster_list, gene2tr){
  
  require(tidyverse)
  
  message(paste("Total no. of clusters:", length(cluster_list), sep = " "))
  message(paste("Total isoforms in clusters:", unlist(cluster_list) %>% length), sep = " ")
  
  # filter out transcripts from genes with only one isoform left
  
      # list of all clustered transcripts
      clustered_tr <- unlist(cluster_list) %>% unname
      # split/group transcripts by gene IDs
      clustered_g <- split(clustered_tr, gene2tr[match(clustered_tr, gene2tr$transcript_id),]$gene_id)
      # count no. of transcripts per gene in the clusters
      ntr <- purrr::map_int(clustered_g, length)
      # find genes with > 1 isoform clustered
      gmulti <- ntr[ntr > 1]
      # find transcripts from genes with > 1 isoform clustered
      multi_tr <- unlist(clustered_g[names(gmulti)])
      # filter clusters to keep only transcripts with multiple same-gene transcripts
      clusters_multi <- purrr::map(cluster_list, ~(.[. %in% multi_tr]))
  
  message(paste("Isoforms clustered after coordination filter:", unlist(clusters_multi) %>% length, sep = " "))
      
  # filter out transcripts from genes with all isoforms in same cluster

      # convert clusters to gene IDs
      clusters_multi.gene <- purrr::map(clusters_multi, ~(gene2tr[match(., gene2tr$transcript_id),]$gene_id))
      
      # find no. of clusters where each gene has isoforms
      gene_distribution <- purrr::map(clusters_multi.gene, ~(names(gmulti) %in% .)) %>% bind_rows %>% t
      colnames(gene_distribution) <- names(gmulti)
      gene_distribution <- colSums(gene_distribution)
      # find differentially spliced genes (i.e. isoforms in more than one cluster)
      genes_ds <- gene_distribution[gene_distribution > 1]
      
      # filter clusters to only keep transcripts from ds genes
      tr_ds.idx <- purrr::map(clusters_multi.gene, ~(which(. %in% names(genes_ds))))
      clusters_ds <- purrr::map2(clusters_multi, tr_ds.idx, ~(.x[.y]))
  
  message(paste("Isoforms clustered after differential splicing filter:", unlist(clusters_ds) %>% length, sep = " "))
  
  return(clusters_ds)
}


##### FUNCTION TO SCALE DATA #####
scale_range <- function(data){
  
  require(tidyverse)
  
  # handle rownames
  id <- str_detect(colnames(data), "transcript")
  
  if(sum(id) == 1){
    
    name_id <- colnames(data[,which(id)])
    data <- column_to_rownames(data, var = name_id)
  }
  
  # scale transript expression
  tr_center <- apply(data, 1, mean)
  
  tr_max <- apply(data, 1, max)
  tr_min <- apply(data, 1, min)
  tr_range <- tr_max - tr_min
  
  data_scaled <- apply(data, 2, function(x) (x - tr_center)/tr_range)
  
  return(data_scaled %>% as.data.frame)
}


##### FUNCTION TO PLOT CLUSTER #####
plot_avg_expr <- function(data, tr_names, cell_types, plot_title = NULL, return = FALSE,
                          labels = NULL){
  
  require(plotrix)
  require(tidyverse)
  require(cowplot)
  
  # handle rownames
  id <- str_detect(colnames(data), "transcript")
  
  if(sum(id) == 1){
    
    name_id <- colnames(data[,which(id)])
    data <- column_to_rownames(data, var = name_id)
  }
  
  # calculate mean by cell type
  splt <- data[tr_names, ] %>% t %>% as.data.frame %>% split(cell_types)
  means <- purrr::map(splt, colMeans) %>% purrr::map(enframe, name = "transcript_id", value = "expression")
  # format as data frame
  means <- bind_rows(means, .id = "cell_type")
  
  # calculate standard error by cell type
  errors <- purrr::map(splt, ~(apply(., 2, std.error) %>% as.data.frame))
  # format as data.frame
  errors <- purrr::map(errors, rownames_to_column) %>% bind_rows
  # add standard error column to means object
  means <- mutate(means, error = errors$.)
  # format cell type factor
  if(is.null(labels) == FALSE){
    means$cell_type <- factor(means$cell_type,
                              levels = unique(means$cell_type) %>% sort(), labels = labels)
  } else if(is.null(labels) == TRUE){
    means$cell_type <- factor(means$cell_type,
                              levels = unique(means$cell_type) %>% sort())
  }
  
  
  # return means or plot
  if(return == TRUE){
    return(means)
    
  } else {
    ggplot(means, aes(x = cell_type, y = expression, colour = transcript_id, group = transcript_id)) + 
      ggtitle(plot_title) +
      geom_line() + geom_point() + 
      geom_errorbar(aes(ymin = expression - error, ymax = expression + error), width = 0.1) + 
      ylab("Mean expression (scaled)") + xlab("Cell type") + 
      theme_set(theme_cowplot()) +
      theme(legend.title = element_blank(), legend.position = "none")
  }
}

##### FUNCTION TO CALCULATE CLUSTER MEAN PROFILE #####
calc_avg_profile <- function(data, tr_names, cell_types, plot_title = NULL, return = FALSE){
  
  require(tidyverse)
  require(cowplot)
  
  # handle rownames
  id <- str_detect(colnames(data), "transcript")
  
  if(sum(id) == 1){
    name_id <- colnames(data[,which(id)])
    data <- column_to_rownames(data, var = name_id)
  }
  
  # all expression values data frame
  clust_avg <- data[tr_names,] %>% t %>% as.data.frame %>% split(cell_types) %>% 
    purrr::map(~(colMeans(.) %>% enframe(name = "transcript_id", value = "ct_mean"))) %>% 
    bind_rows(.id = "cell_type")
  # mean profile
  avg <- aggregate(clust_avg$ct_mean, list(clust_avg$cell_type), mean)
  sd <- aggregate(clust_avg$ct_mean, list(clust_avg$cell_type), sd)

  # format mean profile as data frame
  avg_sil <- tibble(tr = rep("silhouette", length(nrow(clust_avg))), 
                    mean = avg$x, sd = sd$x, cell_type = avg$Group.1)
  # format cell type factor
  #avg_sil$cell_type <- factor(avg_sil$cell_type, 
  #                         levels = unique(avg_sil$cell_type) %>% sort(), 
  #                           labels = c("Astr", "End", "GABA", "Glut", "Micr", "Oligo", "OPC"))
  
  if(return == TRUE){
    return(avg_sil)
  } else {
    # plot all transcripts + mean profile
    avg_all <- plot_avg_expr(data, clust, cell_types, return = TRUE)
    avg <- bind_rows(avg_all, avg_sil)
    avg <- avg %>% dplyr::mutate(category = factor(c(rep("tr", nrow(avg_all)), rep("sil", nrow(avg_sil))),
                                                   levels = c("tr", "sil"),
                                                   labels = c("Transcript", "Mean profile")))
    ggplot(avg, aes(x = cell_type, y = value, colour = category, group = transcript)) + geom_line() + geom_point() + 
      # geom_errorbar(aes(ymin = value - error, ymax = value + error), width = 0.1) + 
      ylab("Mean expression (scaled)") + xlab("Cell type") + theme(legend.title = element_blank()) +
      scale_color_manual(values = c("grey", "black")) + theme_set(theme_cowplot())
  }
  
}




