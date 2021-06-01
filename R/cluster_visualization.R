##### FUNCTION TO SCALE DATA #####

#' @export
scale_range <- function(data){

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

#' @export
plot_avg_expr <- function(data, tr_names, cell_types, plot_title = NULL, return = FALSE,
                          labels = NULL){

  # handle rownames
  id <- str_detect(colnames(data), "transcript")

  if(sum(id) == 1){

    name_id <- colnames(data[,which(id)])
    data <- column_to_rownames(data, var = name_id)
  }

  # calculate mean by cell type
  splt <- data[tr_names, ] %>% t %>% as.data.frame %>% split(cell_types)
  means <- map(splt, colMeans) %>% map(enframe, name = "transcript_id", value = "expression")
  # format as data frame
  means <- bind_rows(means, .id = "cell_type")

  # calculate standard error by cell type
  errors <- map(splt, ~(apply(., 2, plotrix::std.error) %>% as.data.frame))
  # format as data.frame
  errors <- map(errors, rownames_to_column) %>% bind_rows
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
      theme(legend.title = element_blank(), legend.position = "none")
  }
}


##### FUNCTION TO CALCULATE CLUSTER MEAN PROFILE #####

#' @export
calc_avg_profile <- function(data, tr_names, cell_types, plot_title = NULL, return = FALSE){

  # handle rownames
  id <- str_detect(colnames(data), "transcript")

  if(sum(id) == 1){
    name_id <- colnames(data[,which(id)])
    data <- column_to_rownames(data, var = name_id)
  }

  # all expression values data frame
  clust_avg <- data[tr_names,] %>% t %>% as.data.frame %>% split(cell_types) %>%
    map(~(colMeans(.) %>% enframe(name = "transcript_id", value = "ct_mean"))) %>%
    bind_rows(.id = "cell_type")
  # mean profile
  avg <- stats::aggregate(clust_avg$ct_mean, list(clust_avg$cell_type), mean)
  sd <- stats::aggregate(clust_avg$ct_mean, list(clust_avg$cell_type), sd)

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
    avg <- avg %>% mutate(category = factor(c(rep("tr", nrow(avg_all)), rep("sil", nrow(avg_sil))),
                                            levels = c("tr", "sil"),
                                            labels = c("Transcript", "Mean profile")))
    ggplot(avg, aes(x = cell_type, y = value, colour = category, group = transcript)) + geom_line() + geom_point() +
      # geom_errorbar(aes(ymin = value - error, ymax = value + error), width = 0.1) +
      ylab("Mean expression (scaled)") + xlab("Cell type") + theme(legend.title = element_blank()) +
      scale_color_manual(values = c("grey", "black"))
  }

}
