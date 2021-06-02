#' @import purrr
#' @import dplyr
#' @import tibble
#' @import ggplot2



#' @title Scale single-cell expression data for cluster visualization
#'
#' @description This function performs \strong{range scaling} on isoform-level
#' single-cell counts to assist visualization of isoform clusters obtained with
#' \code{acorde}
#'
#' @param data A data.frame or tibble including isoforms as rows and cells
#' as columns. Isoform IDs can be included as row names (data.frame) or as an
#' additional column (tibble).
#'
#' @param isoform_col When a tibble is provided in \code{data}, a character
#' object indicating the name of the column where isoform IDs are specified.
#' Otherwise, isoform identifiers will be assumed to be defined as rownames,
#' and this argument will not need to be provided.
#'
#' @details The purpose of scaling is to be able to jointly visualize isoform
#' expression trends for all members of a cluster, independently of each isoform's
#' absolute expression level.
#'
#' For each isoform, counts are first \strong{centered} by substracting the isoform mean
#' across all cell types. Then, the \strong{expression range} is computed as the difference
#' between the maximum and minimum count values of the isoform. Of note, this
#' range is often equivalent to the maximum counts, since most isoforms show
#' minimum count values of zero. \strong{Scaled} counts are then computed by diving centered
#' counts by the expression range.
#'
#' @return A data.frame containing the scaled counts, with cell IDs as column
#' names and isoform IDs as row names.
#'
#' @export
scale_range <- function(data, isoform_col = NULL){

  # handle rownames
  if(is.null(isoform_col) == FALSE){
    data <- column_to_rownames(data, isoform_col)
  }

  # scale trancsript expression
  tr_center <- apply(data, 1, mean)

  tr_max <- apply(data, 1, max)
  tr_min <- apply(data, 1, min)
  tr_range <- tr_max - tr_min

  data_scaled <- apply(data, 2, function(x) (x - tr_center)/tr_range)

  return(data_scaled %>% as.data.frame)
}



#' @title Calculate isoform-level summary metrics for a cluster
#'
#' @description This function computes cell type mean expression and standard
#' error for all the isoforms in a cluster.
#'
#' @param data A data.frame or tibble including single-cell data, with
#' isoforms as rows and cells as columns. Isoform IDs can be included as row
#' names (data.frame) or as an additional column (tibble). Expression should be
#' previously \strong{scaled} (see \code{\link{scale_range}}).
#'
#' @param isoform_ids A character vector including the IDs of the isoforms included
#' in the cluster.
#'
#' @param id_table A data frame including two columns named \code{cell} and \code{cell_type},
#' providing the correspondence between cell ID and cell type should be
#' provided. The number of rows should be equal to the total number of
#' cell columns in \code{data}, and the order of the \code{cell} column should
#' match column (i.e. cell) order in \code{data}.
#'
#' @param isoform_col When a tibble is provided in \code{data}, a character
#' object indicating the name of the column where isoform IDs are specified.
#' Otherwise, isoform identifiers will be assumed to be defined as rownames,
#' and this argument will not need to be provided.
#'
#' @return A tibble containing four columns: isoform IDs, mean expression values
#' in each of the cell types, cell type names, and standard error of cell type-level
#' expression.
#'
#' @export
calculate_cluster_ctmeans <- function(data, isoform_ids, id_table,
                                      isoform_col = NULL){

  # handle rownames
  if(is.null(isoform_col) == FALSE){
    data <- column_to_rownames(data, isoform_col)
  }

  # calculate mean by cell type
  splt <- data[isoform_ids, ] %>% t %>% as.data.frame %>% split(id_table$cell_type)
  means <- map(splt, colMeans) %>% map(enframe, name = "transcript_id", value = "mean_expr")
  # format as data frame
  means <- bind_rows(means, .id = "cell_type")

  # calculate standard error by cell type
  errors <- map(splt, ~(apply(., 2, plotrix::std.error) %>% as.data.frame))
  # format as data.frame
  errors <- map(errors, rownames_to_column) %>% bind_rows
  # add standard error column to means object
  means <- mutate(means, error = errors$.)

  # return means
  return(means)
}



#' @title Plot cell type expression of all isoforms in a cluster
#'
#' @description This function uses the output of \code{\link{calculate_cluster_ctmeans}}
#' to create a plot including all of the isoforms assigned to a cluster, enabling comparison
#' of their expression pattern across cell types.
#'
#' @param data A tibble including cell type-level mean expression for all isoforms
#' in a cluster (output by \code{\link{calculate_cluster_ctmeans}}).
#'
#' @param plot_title A character object providing a title for the plot.
#'
#' @param ct_labels A character vector including a plot label for each of
#' the cell types defined in \code{data}.
#'
#' @return A plot object generated by \code{\link[ggplot2]{ggplot}}.
#'
#' @export
plot_cluster_ctmeans <- function(data, plot_title = NULL, ct_labels = NULL){

  # format cell type factor
  if(is.null(ct_labels) == FALSE){
    data$cell_type <- factor(data$cell_type,
                              levels = unique(data$cell_type) %>% sort(), labels = ct_labels)
  } else if(is.null(ct_labels) == TRUE){
    data$cell_type <- factor(data$cell_type,
                              levels = unique(data$cell_type) %>% sort())
  }

  # plot
  p <- ggplot(data, aes(x = cell_type, y = mean_expr, colour = transcript_id,
                        group = transcript_id)) +
    ggtitle(plot_title) +
    geom_line() + geom_point() +
    geom_errorbar(aes(ymin = mean_expr - error, ymax = mean_expr + error), width = 0.1) +
    ylab("Mean expression (scaled)") + xlab("Cell type") +
    theme(legend.title = element_blank(), legend.position = "none")

  return(p)
}



#' @title Calculate mean profile of an isoform cluster
#'
#' @description This function computes mean cluster expression across all isoforms,
#' as well as the standard deviation. First, cell type-level mean expression
#' is computed for each isoform. Then, the average value of all cell type means
#' is computed across isoforms.
#'
#' @param data A data.frame or tibble including single-cell data, with
#' isoforms as rows and cells as columns. Isoform IDs can be included as row
#' names (data.frame) or as an additional column (tibble). Expression should be
#' previously \strong{scaled} (see \code{\link{scale_range}}).
#'
#' @param isoform_ids A character vector including the IDs of the isoforms included
#' in the cluster.
#'
#' @param id_table A data frame including two columns named \code{cell} and \code{cell_type},
#' providing the correspondence between cell ID and cell type should be
#' provided. The number of rows should be equal to the total number of
#' cell columns in \code{data}, and the order of the \code{cell} column should
#' match column (i.e. cell) order in \code{data}.
#'
#' @param isoform_col When a tibble is provided in \code{data}, a character
#' object indicating the name of the column where isoform IDs are specified.
#' Otherwise, isoform identifiers will be assumed to be defined as rownames,
#' and this argument will not need to be provided.
#'
#' @return A tibble including three main columns: cell type name, global
#' cluster mean expression in each cell type, and standard deviation computed for
#' the cell type mean expression of all isoforms.
#'
#' @export
calculate_cluster_profile <- function(data, isoform_ids, id_table,
                                      isoform_col = NULL){

  # handle rownames
  if(is.null(isoform_col) == FALSE){
    data <- column_to_rownames(data, isoform_col)
  }

  # all expression values data frame
  clust_avg <- data[isoform_ids,] %>% t %>% as.data.frame %>% split(id_table$cell_type) %>%
    map(~(colMeans(.) %>% enframe(name = "transcript_id", value = "ct_mean"))) %>%
    bind_rows(.id = "cell_type")
  # mean profile
  avg <- stats::aggregate(clust_avg$ct_mean, list(clust_avg$cell_type), mean)
  sd <- stats::aggregate(clust_avg$ct_mean, list(clust_avg$cell_type), stats::sd)

  # format mean profile as data frame
  avg_sil <- tibble(tr = rep("silhouette", length(nrow(clust_avg))),
                    mean = avg$x, sd = sd$x, cell_type = avg$Group.1)

  # return mean profile
  return(avg_sil)
}



#### PLOT CLUSTER MEAN PROFILE ####
#' @title Plot mean profile of an isoform cluster
#'
#' @description This function uses the output of \code{\link{calculate_cluster_profile}}
#' to plot the mean scaled expression of all isoforms in a cluster, enabling a
#' summarized vision of the cluster's global pattern.
#'
#' @param data A tibble including a global mean expression value for each cell type,
#' computed across all isoforms in the cluster (output by
#' \code{\link{calculate_cluster_ctmeans}}).
#'
#' @param plot_title A character object providing a title for the plot.
#'
#' @param ct_labels A character vector including a plot label for each of
#' the cell types defined in \code{data}.
#'
#' @return A plot object generated by \code{\link[ggplot2]{ggplot}}.
#'
#' @export
plot_cluster_profile <- function(data, plot_title = NULL, ct_labels = NULL){

  # format cell type factor
  if(is.null(ct_labels) == FALSE){
    data$cell_type <- factor(data$cell_type,
                             levels = unique(data$cell_type) %>% sort(), labels = ct_labels)
  } else if(is.null(ct_labels) == TRUE){
    data$cell_type <- factor(data$cell_type,
                             levels = unique(data$cell_type) %>% sort())
  }

  # add upper and lower columns for plotting ribbon
  data <- data %>% mutate(upper = mean+sd, lower = mean-sd)

  # plot
  p <- ggplot(data) + ggtitle(plot_title) +
    geom_ribbon(aes(x = cell_type, ymax = upper, ymin = lower, group = 1),
                fill = "gray80") +
    geom_line(aes(x = cell_type, y = mean, group = 1),
              size = 1.5, color = "red3") +
    xlab("Cell type") + ylab("Scaled counts") +
    theme(axis.text.x = element_text(angle = 90))

  return(p)
}

