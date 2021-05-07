##### FUNCTION TO FILTER BY DIFFERENTIAL EXPRESSION OF ISOFORMS #####

# data: SCE object containing counts and cell type labels in a cell_type column

cell_type_DE <- function(data, pvalue_filter = 0.05, mode = c("edgeR", "DESeq2", "both"),
                         compute_weights = TRUE, offset = NULL){
  
  require(tidyverse)
  
  mode <- match.arg(mode)
  
  # calculate cell-level weights for each transcript
  if(compute_weights == TRUE){
    
    require(zinbwave)
    message("Computing cell-level weights using ZINBWaVe...")
    
    data <- zinbwave(data, observationalWeights = TRUE)
    
  } else {
    message("Warning: compute_weights = FALSE. Weights will be required to be pre-computed and stored in the SCE object.")
  }
  
  # DE WITH EDGER (if indicated by user)
  if(mode == "edgeR" || mode == "both"){
    
    require(edgeR)
    message("Running DE across cell types with edgeR...")
    
    # create DGEList object for edgeR and add missing elements for DE
    dge <- DGEList(assay(data))
    if(is.null(offset) == FALSE){
      dge$offset <- offset
    }
    dge <- calcNormFactors(dge)
    design <- model.matrix(~cell_type, data = colData(data))
    dge$weights <- assay(data, "weights")
    
    # glm fit and test
    dge <- estimateDisp(dge, design)
    fit <- glmFit(dge, design)
    lrt <- glmWeightedF(fit, coef = 2:length(unique(data$cell_type)))
    
    # extract and filter results by p-value threshold
    edger_results <- lrt$table %>% rownames_to_column("transcript")
    edger_sig.tr <- filter(edger_results, padjFilter < pvalue_filter) %>% as_tibble()
  }
  
  
  # DE WITH DESEQ2 (if indicated by user)
  if(mode == "DESeq2" || mode == "both"){
    
    require(DESeq2)
    message("Running DE across cell types with DESeq2...")
    
    # create DESeqDataSet object for DESeq2
    dds <- DESeqDataSet(data, design = ~cell_type)
    
    if(is.null(offset) == FALSE){
      normFactors <- exp(-1 * offset)
      normFactors <- normFactors / exp(rowMeans(log(normFactors)))
      normalizationFactors(dds) <- normFactors
    }
    
    # run DESeq() using indications from zinbwave package
    dds <- DESeq(dds, sfType = "poscounts", minReplicatesForReplace = Inf, parallel = TRUE)
    
    # extract and filter results
    deseq_results <- results(dds, independentFiltering = FALSE) %>% as.data.frame
    deseq_sig.tr <- deseq_results %>% rownames_to_column("transcript") %>% as_tibble %>% filter(padj < pvalue_filter) 
  }
  
  
  # OUTPUT
  if(mode == "both"){
    return(list(edgeR = edger_sig.tr, DESeq2 = deseq_sig.tr))
    
  } else if(mode == "edgeR"){
    return(edger_sig.tr)
    
  } else if(mode == "DESeq2"){
    return(deseq_sig.tr)
    
  }
}

##### FUNCTION TO FILTER BY LOW EXPRESSION IN PROPORTION OF A CELL TYPE ####
filter_ctzeros <- function(data, cell_types, ct_proportion = 0.2, rownames = "transcript_id"){
  
  # handle rownames
  data <- data %>% as.data.frame %>% column_to_rownames(rownames)
  
  # test number of zeros in each cell type
  split <- split(data %>% t %>% as.data.frame, cell_types)
  test_zero <- map(split, ~(. > 0))
  
  # compare to cell type proportion threshold
  lgl <- map(test_zero, ~(colSums(.) >= nrow(.)*ct_proportion))
  allct_lgl <- bind_rows(lgl)
  
  # no. of zeros should be higher than proportion in at least one cell type
  final_lgl <- rowSums(allct_lgl) >= 1
  
  return(final_lgl)
}


# FUNCTION TO RUN ONE DOWNSAMPLING ITERATION
run_downsampling <- function(data, cell_types, cell_no = 45){
  
  message("Downsampling data...")
  
  require(zinbwave)
  
  # randomly sample 50 cells from each neural type
  gaba_ids <- filter(cell_types, cell_type == "GABA-ergic Neuron") %>% 
    select(run) %>% unlist %>% sample(45)
  glut_ids <- filter(cell_types, cell_type == "Glutamatergic Neuron") %>% 
    select(run) %>% unlist %>% sample(45)
  noneur_ids <- filter(cell_types, cell_type != "Glutamatergic Neuron" & 
                         cell_type != "GABA-ergic Neuron") %>% 
    select(run) %>% unlist
  
  # subset expression matrix
  data_down <- data[,c("transcript_id", gaba_ids, glut_ids, noneur_ids)]
  # subset metadata and reorder columns in matrix
  cell_types <- filter(cell_types, run %in% colnames(data_down %>% select(-transcript_id)))
  data_down <- data_down[, c("transcript_id", cell_types$run)]
  
  # repeat filtering: ensure high proportion of expression in at least one cell type
  keep_tr <- filter_ctzeros(data_down, cell_types$cell_type, ct_proportion = 0.25)
  data_down <- data_down[keep_tr,]
  message(paste0("Total genes to be tested for DE: ", nrow(data_down)))
  
  # create SCE object
  sce <- SingleCellExperiment(assays = list(counts = column_to_rownames(data_down, "transcript_id") %>% 
                                              as.matrix %>% round),
                              colData = cell_types)
  
  # calculate weights
  message("Calculating ZINBWaVE weights...")
  sce <- zinbwave(sce, observationalWeights = TRUE, BPPARAM = MulticoreParam(6))
  
  return(sce)
}