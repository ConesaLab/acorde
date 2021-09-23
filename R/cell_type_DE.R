#' @import dplyr
#' @import purrr
#' @import tibble
#' @importFrom Rdpack reprompt



#' @title Test differential expression of isoforms across cell types
#'
#' @description This function provides a wrapper to compute Differential
#' Expression across multiple cell types on a single-cell dataset, using
#' standard bulk tools edgeR and DESeq2 and ZINBWaVE weights, as described
#' by Van den Berge et al. (see references below).
#'
#' @param data A \code{\link[SingleCellExperiment]{SingleCellExperiment}} (SCE)
#' object including isoform-level counts on the \code{assay()} slot and a
#' \code{cell_type} column indicating cell type correspondence in the
#' \code{colData()} slot.
#' @param mode A character value indicating whether \code{"edgeR"},
#' \code{"DESeq2"} or \code{"both"} DE methods are to be run.
#' @param compute_weights A logical value. If \code{TRUE}, cell-level weights
#' will be computed using ZINBWaVE (see \code{\link[zinbwave]{zinbwave}}).
#' Alternatively, weights may be pre-computed and stored in the \code{weights()}
#' slot in the SCE object.
#' @param AdjPvalue A numeric value. If provided, filtering of DE results by
#' adjusted p-value is enabled.
#' @param maxFC A numeric value. If provided, filtering of DE results by
#' maximum fold-change (maxFC) across cell types is enabled. Note that FC filter
#' is based on edgeR's output (FC between all conditions), and therefore this
#' option requires modes \code{"edgeR"} or \code{"both"} to be enabled.
#' @param offset Optionally, a matrix containing offset information to be
#' considered for bias correction in edgeR and/or DESeq2 DE tests.
#'
#' @return If \code{mode = "both"}, a named list of length 2 containing the results
#' (filtered or not, depending on whether \code{AdjPvalue} and \code{maxFC} are
#' provided) of testing Differential Expression using edgeR and DESeq2. If
#' \code{maxFC} is provided, columns \code{max_FC} and \code{max_log2FC} are added
#' to the results tables. If the \code{mode} argument was used to select only
#' one method, the function returns a data.frame object containing the results
#' table obtained by that method.
#'
#' @references
#' \insertRef{Robinson2009}{acorde}
#'
#' \insertRef{Love2014}{acorde}
#'
#' \insertRef{VandenBerge2018}{acorde}
#'
#' @export
cell_type_DE <- function(data, mode = c("edgeR", "DESeq2", "both"),
                         compute_weights = TRUE,
                         AdjPvalue = NULL,
                         maxFC = NULL,
                         offset = NULL){

  mode <- match.arg(mode)

  # check install
  if(!requireNamespace("zinbwave", quietly = TRUE)){
    stop("package 'zinbwave' is required but not installed.")
  }
  if(!requireNamespace("DESeq2", quietly = TRUE)){
    stop("package 'DESeq2' is required but not installed.")
  }
  if(!requireNamespace("edgeR", quietly = TRUE)){
    stop("package 'edgeR' is required but not installed.")
  }
  if(!requireNamespace("SingleCellExperiment", quietly = TRUE)){
    stop("package 'SingleCellExperiment' is required but not installed.")
  }

  # calculate cell-level weights for each transcript
  if(compute_weights == TRUE){

    message("Computing cell-level weights using ZINBWaVe...")

    data <- zinbwave::zinbwave(data, observationalWeights = TRUE,
                               BPPARAM = BiocParallel::MulticoreParam(6))

  } else {
    message("Note: compute_weights = FALSE.
            Weights will be required to be pre-computed and stored in the
            SCE object.")
  }

  # DE WITH EDGER (if indicated by user)
  if(mode == "edgeR" || mode == "both"){

    message("Running DE across cell types with edgeR...")

    # create DGEList object for edgeR and add missing elements for DE
    dge <- edgeR::DGEList(SummarizedExperiment::assay(data))
    if(is.null(offset) == FALSE){
      dge$offset <- offset
    }
    dge <- edgeR::calcNormFactors(dge)
    design <- stats::model.matrix(~cell_type, data = SummarizedExperiment::colData(data))
    dge$weights <- SummarizedExperiment::assay(data, "weights")

    # glm fit and test
    dge <- edgeR::estimateDisp(dge, design)
    fit <- edgeR::glmFit(dge, design)
    lrt <- zinbwave::glmWeightedF(fit, coef = 2:length(unique(data$cell_type)))

    # extract results
    edger_results <- lrt$table %>% rownames_to_column("transcript")

    # if provided, filter by p-value threshold
    if(is.null(AdjPvalue) == FALSE){

      message("Applying AdjPvalue filter...")

      edger_sig.tr <- dplyr::filter(edger_results,
                                    padjFilter < AdjPvalue) %>%
        tibble::as_tibble()

    }else{

      message("Note: p-value threshold not provided, full edgeR results will be returned.
              To enable filtering, please set AdjPvalue")
      edger_sig.tr <- edger_results

    }

    # if provided, filter by FC threshold
    if(is.null(maxFC) == FALSE){

      message("Filtering DE transcripts by maximum fold-change between cell types...")

      # compute maximum FC and log2FC
      fc_edger <- edger_sig.tr %>%
        dplyr::select(tidyselect::starts_with("logFC")) %>%
        dplyr::rowwise() %>%
        abs() %>%
        dplyr::mutate(max_log2FC = max(dplyr::across(tidyselect::everything())),
                      max_FC = 2^max_log2FC) %>%
        dplyr::select(max_log2FC, max_FC)

      # add fc to results
      edger_sig.tr <- dplyr::bind_cols(edger_sig.tr, fc_edger)

      # filter
      edger_sig.tr <- dplyr::filter(edger_sig.tr,
                                    max_FC > maxFC)

    }else{
      message("Note: maximum fold-change threshold not provided, results will not
              be filtered by cell-type fold-change. To enable filtering, please set maxFC.")
    }
  }


  # DE WITH DESEQ2 (if indicated by user)
  if(mode == "DESeq2" || mode == "both"){

    message("Running DE across cell types with DESeq2...")

    # create DESeqDataSet object for DESeq2
    dds <- DESeq2::DESeqDataSet(data, design = ~cell_type)

    if(is.null(offset) == FALSE){
      normFactors <- exp(-1 * offset)
      normFactors <- normFactors / exp(rowMeans(log(normFactors)))
      DESeq2::normalizationFactors(dds) <- normFactors
    }

    # run DESeq() using indications from zinbwave package
    dds <- DESeq2::DESeq(dds, sfType = "poscounts", minReplicatesForReplace = Inf, parallel = TRUE)

    # extract results
    deseq_results <- DESeq2::results(dds, independentFiltering = FALSE) %>% as.data.frame

    # filter by p-value threshold
    if(is.null(AdjPvalue) == FALSE){

      message("Applying AdjPvalue filter...")

      deseq_sig.tr <- deseq_results %>% rownames_to_column("transcript") %>% as_tibble %>%
        filter(padj < AdjPvalue)

    }else{

      message("Note: p-value threshold not provided, full DESeq2 results will be returned.
              To enable filtering, please set AdjPvalue")
      deseq_sig.tr <- deseq_results

    }

    # filter by FC threshold
    if(is.null(maxFC) == FALSE && mode == "both"){

      # compute maximum FC and log2FC
      fc_deseq <- deseq_sig.tr %>%
        dplyr::select(starts_with("logFC")) %>%
        dplyr::rowwise() %>%
        abs() %>%
        dplyr::mutate(max_log2FC = max(dplyr::across(tidyselect::everything())),
                      max_FC = 2^max_log2FC) %>%
        dplyr::select(max_log2FC, max_FC)

      # add fc to results
      deseq_sig.tr <- dplyr::bind_cols(deseq_sig.tr, fc_deseq)

      # filter
      deseq_sig.tr <- dplyr::filter(deseq_sig.tr,
                                    max_FC > maxFC)

    }else{
      message("Note: maximum fold-change threshold not provided, results will not
              be filtered by cell-type fold-change. To enable filtering, please set maxFC.")
    }
  }

  message(paste0("Finished running DE analysis with mode = ", mode, "."))

  # OUTPUT
  if(mode == "both"){
    return(list(edgeR = edger_sig.tr, DESeq2 = deseq_sig.tr))

  } else if(mode == "edgeR"){
    return(edger_sig.tr)

  } else if(mode == "DESeq2"){
    return(deseq_sig.tr)

  }
}



#' @title Filter isoforms by maximum cell type-level proportion of zeros
#'
#' @description For convenience, \code{acorde} includes a function to assist
#' lenient filtering of isoforms based on the proportion of zero expression values
#' across cell types. Based on this criteria, a minimum number of cells must have
#' non-zero expression in at least one cell type.
#'
#' @param data A data.frame or tibble object including isoforms as rows and cells as columns.
#' Isoform IDs can be included as row names (data.frame) or as an additional column (tibble).
#'
#' @param id_table A data frame including two columns named \code{cell} and \code{cell_type},
#' in which correspondence between cell ID and cell type should be
#' provided. The number of rows should be equal to the total number of
#' cell columns in \code{data}, and the order of the \code{cell} column should
#' match column (i.e. cell) order in \code{data}.
#'
#' @param ct_proportion A numeric indicating the minimum proportion of cells with
#' non-zero expression that will be allowed per cell type. Isoforms with a non-zero
#' value proportion above the threshold in at least one cell type will be flagged
#' to be preserved. Defaults to 0.2 (i.e. 20\%).
#'
#' @param isoform_col When a tibble is provided in \code{data}, a character object
#' indicating the name of the column in which isoform IDs are specified.
#' Otherwise, isoform identifiers will be assumed to be defined as rownames,
#' and this argument will not need to be provided.
#'
#' @return A logical vector including one entry per isoform in \code{data}. Isoforms
#' meeting the sparsity criteria will have a value of \code{TRUE}, and otherwise
#' be labeled as \code{FALSE}. This logical vector can then be used to filter
#' isoforms, i.e. the rows in \code{data}.
#'
#' @export
detect_sparse <- function(data, id_table, ct_proportion = 0.2, isoform_col = NULL){

  # handle rownames
  if(is.null(rownames) == FALSE){
    data <- data %>% as.data.frame %>% column_to_rownames(isoform_col)
  }

  # test number of zeros in each cell type
  split <- split(data %>% t %>% as.data.frame, id_table$cell_type)
  test_zero <- map(split, ~(. > 0))

  # compare to cell type proportion threshold
  lgl <- map(test_zero, ~(colSums(.) >= nrow(.)*ct_proportion))
  allct_lgl <- bind_rows(lgl)

  # no. of non-zero values should be higher than proportion in at least one cell type
  final_lgl <- rowSums(allct_lgl) >= 1

  return(final_lgl)
}


# FUNCTION TO RUN ONE DOWNSAMPLING ITERATION

#' @title Run one downsampling iteration on Tasic et al. dataset
#'
#' @description Function designed to pack random cell downsampling code to ease
#' running multiple iterations of this process during \code{acorde} benchmarking.
#' Each of these iterations includes three steps: randomly selecting a number of
#' cells form the specified cell types, subsetting the expression matrix and
#' creating an SCE object including ZINBWaVE weights for DE testing
#' (see \code{\link{cell_type_DE}}).
#'
#' @param data A data.frame or tibble object including isoforms as rows and cells
#' as columns. Isoform IDs should be included in an independent column, not defined
#' as \code{rownames}.
#' @param id_table A data frame including two columns named \code{cell} and
#' \code{cell_type}, in which correspondence between cell ID and cell type should be
#' provided.
#' @param downsampling_ct A character vector including one or more cell type names
#' (matching those in \code{id_table$cell_type}) to be targeted by downsampling.
#' @param cell_no A numeric indicating the number of cells to be randomly sampled
#' during downsampling. Should be the same for all targeted cell types.
#' @param isoform_col Name of the column in \code{data} that contains isoform IDs.
#' Otherwise, isoform identifiers will be assumed to be defined as rownames,
#' and this argument will not need to be provided.
#'
#' @return An SCE object containing the data after downsampling in
#' \code{assay(counts = data)} and \code{id_table} as
#' metadata in the \code{colData()} slot. This format corresponds to the input of
#' \code{\link{cell_type_DE}}.
#'
#' @export
run_downsampling <- function(data, id_table, downsampling_ct,
                             cell_no, isoform_col = NULL){

  message("Downsampling data...")

  # check install
  if(!requireNamespace("zinbwave", quietly = TRUE)){
    stop("package 'zinbwave' is required but not installed.")
  }
  if(!requireNamespace("SingleCellExperiment", quietly = TRUE)){
    stop("package 'SingleCellExperiment' is required but not installed.")
  }

  # handle rownames and data type
  if(is.null(isoform_col) == TRUE){
    data <- data %>% as.data.frame %>% rownames_to_column("transcript")
    isoform_col <- "transcript"
  }

  # randomly sample 50 cells from each neural type
  down_ids <- map(downsampling_ct,
                  ~filter(id_table, cell_type == .) %>%
                    select(cell) %>%
                    unlist %>% sample(cell_no))
  full_ids <- filter(id_table, !(cell_type %in% downsampling_ct)) %>%
    select(cell_id_column) %>% unlist

  # subset expression data frame
  data_down <- data[,c(isoform_col, unlist(down_ids), full_ids)]
  # subset metadata and reorder columns in data frame
  id_table <- filter(id_table, cell %in% colnames(data_down %>% select(-transcript_id)))
  data_down <- data_down[, c(isoform_col, id_table$cell)]

  # repeat filtering: ensure high proportion of expression in at least one cell type
  keep_tr <- detect_sparse(data_down, id_table$cell_type, ct_proportion = 0.25)
  data_down <- data_down[keep_tr,]
  message(paste0("Total genes to be tested for DE: ", nrow(data_down)))

  # create SCE object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = column_to_rownames(data_down, isoform_col) %>% as.matrix %>% round),
    colData = id_table)

  # calculate weights
  message("Calculating ZINBWaVE weights...")
  sce <- zinbwave::zinbwave(sce, observationalWeights = TRUE,
                            BPPARAM = BiocParallel::MulticoreParam(6))

  return(sce)
}
