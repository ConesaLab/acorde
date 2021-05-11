#' @import purrr
#' @import furrr
#' @import dplyr
#' @import tibble

##### FUNCTION TO FIND SHARED GENES #####
# gene2tr: data frame containing transcript_id and gene_id columns,
#indicating gene-isoform correspondence.

#' @export
find_shared_genes <- function(cluster_list, gene2tr){

  # convert clusters to gene IDs
  cluster_list.gene <- map(cluster_list, ~gene2tr[match(., gene2tr$transcript_id),]$gene_id)

  # find intersections with modified UpSetR fromList() function
  gene_occurrence <- fromList(cluster_list.gene)

  # iterate all combinations of genes and find coincidences across clusters
  gcomb <- utils::combn(rownames(gene_occurrence), 2)
  colnames(gcomb) <- paste0("pair", seq(1, ncol(gcomb)))

  message("Finding co-spliced gene pairs across clusters...")
  check <- future_map_lgl(as_tibble(gcomb), check_gene_pair, gene_occurrence, .progress = TRUE)

  # keep pairs with coincidence
  pairs <- gcomb[,check]

  return(pairs)
}

##### FUNCTION TO CHECK CO-OCCURRENCE OF A GENE PAIR ACROSS CLUSTERS ####
# function called internally by find_shared_genes() to iteratively
# check each gene pair for co-occurrence
check_gene_pair <- function(pair, intersections){

  check <- sum(colSums(intersections[pair, ]) == 2) >= 2
  return(check)
}


##### FUNCTION TO STATISTICALLY TEST ALL SHARED GENES #####

#' @export
test_shared_genes <- function(data, cluster_list, shared_genes, gene2tr, cell_types){

  pair_seq <- seq(1, ncol(shared_genes))

  # run ANOVA test for each gene pair with model fit and test error handling
  pvalues <- future_map(pair_seq,
                        ~tryCatch(make_test(shared_genes[,.],
                                     data, cluster_list, gene2tr, cell_types),
                                  error = function(c){
                                  msg <- conditionMessage(c)
                                  message(paste("Couldn't test pair",
                                          paste0(shared_genes[,.][1], " ", shared_genes[,.][2], ":"),
                                          "glm() error", msg, sep = " "))
                                  NA
                           })
                 )
  return(pvalues)
}

##### FUNCTION TO TEST INTERACTIONS FOR A GENE PAIR ####
# called internally by test_shared_genes to test cluster-dependent
# expression for a pair of shared genes
make_test <- function(pair, data, cluster_list, gene2tr, cell_types){

  # make long matrix with factors and expression
  design <- make_design(data, cluster_list, gene2tr, cell_types, pair)

  # fit glm for all double interactions
  fit <- stats::glm(expression ~ (gene + cluster + cell_type)^2,
             data = design, family = MASS::negative.binomial(theta = 10), maxit = 200)


  # perform test
  adev <- car::Anova(fit, type = 2, contrasts=list(topic=stats::contr.sum, sys=stats::contr.sum))
  # contrasts = list(topic = contr.sum, sys = contr.sum)
  # select p-values
  pvalues <- tibble(adev["cluster:cell_type", "Pr(>Chisq)"], adev["gene:cell_type", "Pr(>Chisq)"])
  colnames(pvalues) <- c("cluster:cell_type", "gene:cell_type")

  return(pvalues)
}


##### FUNCTION TO CREATE DESIGN MATRIX FOR A GENE PAIR ####
# called internally by test_shared_genes to create design matrix for glm for a given pair of shared genes
# pair is a character vector with gene names for the selected pair
make_design <- function(data, cluster_list, gene2tr, cell_types, pair){

  # remake gene_occurrence list
  # convert clusters to gene IDs
  cluster_list.gene <- map(cluster_list, ~gene2tr[match(., gene2tr$transcript_id),]$gene_id)
  # find intersections with modified UpSetR fromList() function
  gene_occurrence <- fromList(cluster_list.gene)

  # find clusters where pair has isoforms
  check <- which(colSums(gene_occurrence[pair, ]) == 2) %>% unname
  # select the clusters (with gene names)
  clust_selec <- cluster_list.gene[check]
  # find transcripts of the genes in each cluster
  idx_list <- map(clust_selec, ~(match(pair, .)))
  tr <- map2(check, idx_list, ~(cluster_list[[.x]][.y])) %>% unlist %>% unname

  # create factor df
      # make data frame with expression of gene pair transcripts
      fouriso_expr <- data %>% filter(transcript_id %in% tr)
      fouriso_expr <- fouriso_expr[match(tr, fouriso_expr$transcript_id),]
      # add gene factor
      gene_factor <- factor(rep(seq(nrow(fouriso_expr)/2), 2))
      fouriso_expr <- mutate(fouriso_expr, gene = gene_factor)
      # add cluster factor
      cluster_factor <- factor(rep(seq(nrow(fouriso_expr)/2), each = 2))
      fouriso_expr <- mutate(fouriso_expr, cluster = cluster_factor)
      # long formatting
      fouriso_expr_long <- tidyr::gather(fouriso_expr, cell_id, expression, -transcript_id, -gene, -cluster)
      fouriso_expr_long <- mutate(fouriso_expr_long, cell_type = rep(cell_types, each = nrow(fouriso_expr))) %>%
        mutate_at("expression", as.integer)

  return(fouriso_expr_long)
}

##### FUNCTION TO CREATE EXPRESSION MATRIX FOR A GENE PAIR ####
# called internally by test_shared_genes to create one-line expression matrix
# to test ct-cluster and ct-gene interactions via a glm approach
make_matrix <- function(data, design){

  # create wide dataession matrix with sample IDs combining cell and transcript
  data <- data %>% filter(transcript_id %in% unique(design$transcript_id))
  data_long <- data %>% tidyr::gather(cell_id, expression, -transcript_id)
  data_wide <- data_long %>% tidyr::pivot_wider(names_from = c(transcript_id, cell_id), values_from = expression)

  # create sample names combining transcript and cell IDs
  design <- design %>% mutate(sample = paste(transcript_id, cell_id, sep = "_")) %>%
    select(sample, gene, cluster, cell_type)

  # reorder columns in wide dataession matrix to match design matrix
  data_wide <- data_wide[,design$sample]

  return(data_wide)
}



