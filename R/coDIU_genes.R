#' @import purrr
#' @import furrr
#' @import dplyr
#' @import tibble


#' @title Detect shared genes for co-Differential Isoform Usage analysis
#'
#' @description This function reports pairs of genes that present co-expressed
#' isoforms given a list of previously-detected isoform clusters.
#' The aim is to enable co-Differential Isoform Usage (co-DIU) analysis on the
#' returned gene pairs, which will be candidates for co-DIU.
#'
#' @param cluster_list A list of character vectors containing isoform IDs.
#' Each element of the list represents a cluster of isoforms.
#' @param gene_tr_table A data.frame or tibble object containing two columns
#' named \code{transcript_id} and \code{gene_id}, indicating gene-isoform
#' correspondence.
#' @param parallel A logical. When \code{TRUE}, parallelization is enabled.
#' The \code{\link{future_map_lgl}} function in the \code{furrr} is used.
#' @param t An integer indicating the number of threads to be used for
#' parallelization. This will be passed to the \code{\link{plan}} function from
#' the \code{future} package via the \code{workers} argument.
#'
#' @details We define coordinated splicing patterns as a situation where
#' post-transcriptional regulation, defined by isoform expression,
#' can be detected independently of transcriptional regulation,
#' i.e. gene-level expression. To detect splicing coordination, we defined
#' co-Differential Isoform Usage (co-DIU) as a pattern where a group of genes
#' shows co-expression of their isoforms, but no co-expression can be detected
#' when only gene expression is considered. In the context of our pipeline,
#' a set of potentially co-DIU genes will have at least two of their isoforms
#' assigned to the same clusters, therefore showing detectable isoform-level
#' co-expression, and suggesting coordinated splicing regulation in that group
#' of genes (see Arzalluz-Luque et al. 2021).
#'
#' @return A matrix containing two rows and as many columns as potentially
#' co-DIU gene pairs detected, that is, genes co-expressing isoforms across
#' two or more clusters. Genes will be represented by the IDs provided in
#' \code{{gene_tr_table}}.
#'
#' @references
#'
#' \insertRef{Arzalluz-Luque2021}{acorde}
#'
#' @export
find_shared_genes <- function(cluster_list, gene_tr_table,
                              parallel = TRUE, t = 4){

  # convert clusters to gene IDs
  cluster_list.gene <- map(cluster_list,
                           ~gene_tr_table[match(., gene_tr_table$transcript),]$gene)

  # find intersections with modified UpSetR fromList() function
  gene_occurrence <- fromList(cluster_list.gene)

  # iterate all combinations of genes and find coincidences across clusters
  gcomb <- utils::combn(rownames(gene_occurrence), 2)
  colnames(gcomb) <- paste0("pair", seq(1, ncol(gcomb)))

  message("Finding co-spliced gene pairs across clusters...")

  if(parallel == TRUE){
    future::plan("multisession", workers = t)
  }
  check <- future_map_lgl(as_tibble(gcomb), check_gene_pair,
                          gene_occurrence)

  # keep pairs with coincidence
  pairs <- gcomb[,check]

  return(pairs)
}


#' @title Check isoform co-expression patterns for a gene pair across clusters
#'
#' @description This function is called internally by \code{\link{find_shared_genes}}
#' to iteratively check each gene pair for isoform co-expression.
#'
#' @param pair A character vector of length 2 containing gene IDs for a candidate
#' gene pair.
#' @param intersections A data.frame including binary presence/absence patterns
#' for genes across clusters. Typically, the output of \code{\link{fromList}}.
#'
#' @return A logical indicating whether the pair of genes under evaluation presents
#' isoform co-expression across any of the analyzed clusters.
check_gene_pair <- function(pair, intersections){

  check <- sum(colSums(intersections[pair, ]) == 2) >= 2
  return(check)
}


#' @title Statistical testing of candidate co-DIU genes
#'
#' @description Pairwise statistical testing of co-Differential Isoform Usage
#' relationships. For a selected set of gene pairs showing co-expression of
#' isoforms across clusters (see \code{\link{find_shared_genes}}), this
#' function tests the significance of the detected co-DIU patterns.
#'
#' @param data A data.frame or tibble object including isoforms as rows and cells
#' as columns. Isoform IDs can be included as row names (data.frame) or as an
#' additional column (tibble).
#'
#' @param cluster_list A list of character vectors containing isoform IDs.
#' Each element of the list represents a cluster of isoforms.
#'
#' @param shared_genes A two-row matrix containing \emph{n} candidate co-DIU
#' gene pairs as column. Typically the result of running
#' \code{\link{find_shared_genes}}.
#'
#' @param gene_tr_table A data.frame or tibble object containing two columns
#' named \code{transcript_id} and \code{gene_id}, indicating gene-isoform
#' correspondence.
#'
#' @param id_table  A data frame including two columns named \code{cell}
#' and \code{cell_type}, in which correspondence between cell ID and cell type
#' should be provided. The number of rows should be equal to the total number of
#' cell columns in \code{data}, and the order of the \code{cell} column should
#' match column (i.e. cell) order in \code{data}.
#'
#' @param isoform_col When a tibble is provided in \code{data}, a character
#' object indicating the name of the column where isoform IDs are specified.
#' Otherwise, isoform identifiers will be assumed to be defined as rownames,
#' and this argument will not need to be provided.
#'
#' @param parallel A logical. When \code{TRUE}, parallelization is enabled.
#' The \code{\link{future_map_lgl}} function in the \code{furrr} is used.
#'
#' @param t An integer indicating the number of threads to be used for
#' parallelization. This will be passed to the \code{\link{plan}} function from
#' the \code{future} package via the \code{workers} argument.
#'
#' @details A set of \strong{potentially co-DIU genes} will have at least two of their
#' isoforms assigned to the same clusters, i.e. show detectable
#' isoform-level co-expression. However, since clustering allows isoforms
#' with slightly variable expression patterns to be clustered together,
#' some isoforms might be assigned to clusters that do not faithfully represent their
#' expression profile, leading to inaccuracies in co-DIU detection.
#' To \strong{avoid false-positive co-DIU genes}, the present function
#' applies a regression model and a statistical
#' test to each of the candidate pair of genes (hereby named
#' gene 1 and gene 2), where at least two of the isoforms of each gene must
#' belong to the same two clusters (hereby named cluster 1 and cluster 2).
#'
#' Briefly, we need to assess whether expression values for the isoforms follow
#' a correct co-DIU pattern, that is, the average profile across cell types of
#' the two isoforms in cluster 1 must be significantly different to the average
#' profile of the two isoforms in cluster 2, indicating distinct expression
#' profiles for the two isoforms of each gene. In addition, the average profile
#' of the two isoforms of gene 1 must not be different to the average profile of
#' the two isoforms of gene 2, indicating that co-expression is only detectable
#' when isoform-level expression is considered.
#'
#' Internally, the function fits a \strong{generalized linear regression model} (GLM) via
#' the \code{\link[stats]{glm}} function, using the \code{\link[MASS]{negative.binomial}}
#' function in the \code{MASS} package to set the error distribution and link
#' function of the model via the \code{family} argument. To test the
#' significance of the \emph{cluster*cell type}
#' and \emph{gene*cell type} interactions (as described above), we calculated
#' \strong{type-II analysis-of-variance} (ANOVA) tables for the model using a
#' \strong{likelihood-ratio chi-square test} using the \code{\link[car]{Anova}} function in
#' the \code{car} package (given the unbalanced design).
#'
#' @return A list containing one \code{tibble} per tested gene pair, as generated
#' by \code{\link{make_test}}. Each tibble will include two columns,
#' \code{cluster:cell_type} and \code{gene:cell_type}, containing the p-value
#' obtained when testing each of these interactions in the type-II
#' ANOVA test.
#'
#' \strong{NOTE:} In some cases the assumptions required for fitting the GLM are not met,
#' and an \code{NA} value is returned instead. These are output to allow users
#' to control for untested gene pairs, but can easily be removed from the output.
#'
#' @seealso For details, see internal functions:
#' \code{\link{make_design}}, \code{\link{make_test}}.
#'
#' @references
#'
#' \insertRef{Venables2002}{acorde}
#'
#' \insertRef{Fox2019}{acorde}
#'
#' @export
test_shared_genes <- function(data, cluster_list, shared_genes, gene_tr_table,
                              id_table, isoform_col = NULL,
                              parallel = TRUE, t = 4){

  # handle rownames and data type
  if(is.null(isoform_col) == TRUE){
    data <- data %>% as.data.frame %>% rownames_to_column("transcript")
  }

  # create vector to iterate
  pair_seq <- seq(1, ncol(shared_genes))

  # run ANOVA test for each gene pair with model fit and test error handling
  if(parallel == TRUE){
    future::plan(multisession, workers = t)
  }

  pvalues <- future_map(pair_seq,
                        ~tryCatch(make_test(shared_genes[,.],
                                     data, cluster_list, gene_tr_table, id_table$cell_type),
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




#' @title Perform co-DIU statistical test for one gene pair
#'
#' @description This function is called iteratively by \code{\link{test_shared_genes}}
#' to test cluster-dependent expression for a single pair of candidate co-DIU genes.
#'
#' @param pair A character vector of length 2 indicating the gene IDs of
#' a gene pair to test.
#'
#' @param data A data.frame or tibble including isoforms as rows and cells
#' as columns. Isoform IDs can be included as row names (data.frame) or as an
#' additional column (tibble).
#'
#' @param cluster_list A list of character vectors containing isoform IDs.
#' Each element of the list represents a cluster of isoforms.
#'
#' @param gene_tr_table A data.frame or tibble containing two columns
#' named \code{transcript} and \code{gene}, indicating gene-isoform
#' correspondence.
#'
#' @param cell_types A character vector including cell type assignments for
#' the cell IDs in \code{data}.
#'
#' @details Arguments are passed from \code{\link{test_shared_genes}}.
#'
#' @return A tibble containing two columns, \code{cluster:cell_type} and
#' \code{gene:cell_type}, containing the p-value obtained in the likelihood ratio
#' chi-square test (see \code{\link{test_shared_genes}}) for each of these two
#' interactions.
#'
#' @references
#'
#' \insertRef{Venables2002}{acorde}
#'
#' \insertRef{Fox2019}{acorde}
make_test <- function(pair, data, cluster_list, gene_tr_table, cell_types){

  # make long matrix with factors and expression
  design <- make_design(pair, data, cluster_list, gene_tr_table, cell_types)

  # fit glm for all double interactions
  fit <- stats::glm(expression ~ (gene + cluster + cell_type)^2,
             data = design, family = MASS::negative.binomial(theta = 10), maxit = 200)

  # perform test
  adev <- car::Anova(fit, type = 2, contrasts=list(topic=stats::contr.sum, sys=stats::contr.sum))
  # select p-values
  pvalues <- tibble(adev["cluster:cell_type", "Pr(>Chisq)"], adev["gene:cell_type", "Pr(>Chisq)"])
  colnames(pvalues) <- c("cluster:cell_type", "gene:cell_type")

  return(pvalues)
}



#' @title Create design matrix for co-DIU test of a pair of genes
#'
#' @description This function is used by \code{\link{test_shared_genes}} to
#' create a design matrix for GLM fitting using information from a pair of co-DIU
#' candidategenes, i.e. expression values of the gene's co-expressed isoforms,
#' cell type labels for cells and cluster labels for isoforms.
#'
#' @param pair A character vector of length 2 indicating the gene IDs of
#' a gene pair to test.
#'
#' @param data A data.frame or tibble including isoforms as rows and cells
#' as columns. Isoform IDs can be included as row names (data.frame) or as an
#' additional column (tibble).
#'
#' @param cluster_list A list of character vectors containing isoform IDs.
#' Each element of the list represents a cluster of isoforms.
#'
#' @param gene_tr_table A data.frame or tibble containing two columns
#' named \code{transcript} and \code{gene}, indicating gene-isoform
#' correspondence.
#'
#' @param cell_types A character vector including cell type assignments for
#' the cell IDs in \code{data}.
#'
#' @details Arguments are passed from \code{\link{test_shared_genes}}.
#'
#' @return A tibble, containing a long-form table with the required factors for
#' GLM fitting and statistical testing.
make_design <- function(pair, data, cluster_list, gene_tr_table, cell_types){

  # remake gene_occurrence list
  # convert clusters to gene IDs
  cluster_list.gene <- map(cluster_list, ~gene_tr_table[match(., gene_tr_table$transcript),]$gene)
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
      fouriso_expr <- data %>% filter(transcript %in% tr)
      fouriso_expr <- fouriso_expr[match(tr, fouriso_expr$transcript),]
      # add gene factor
      gene_factor <- factor(rep(seq(nrow(fouriso_expr)/2), 2))
      fouriso_expr <- mutate(fouriso_expr, gene = gene_factor)
      # add cluster factor
      cluster_factor <- factor(rep(seq(nrow(fouriso_expr)/2), each = 2))
      fouriso_expr <- mutate(fouriso_expr, cluster = cluster_factor)
      # long formatting
      fouriso_expr_long <- tidyr::gather(fouriso_expr, cell_id, expression,
                                         -transcript, -gene, -cluster)
      fouriso_expr_long <- mutate(fouriso_expr_long,
                                  cell_type = rep(cell_types,
                                                  each = nrow(fouriso_expr))) %>%
        mutate_at("expression", as.integer)

  return(fouriso_expr_long)
}
