### Helper functions for ST data analysis


#'
#'
#'
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}


#' merge dataframe with locations
#' internal function for plotting data
#' tp should be a dataframe with pxl_col, pxl_row
merge_with_locations <- function(.data, tp){
  .data <- .data %>%
    merge(tp[c('pxl_col', 'pxl_row')], by = 0) %>%
    tibble::column_to_rownames("Row.names") %>%
#    dplyr::select(-in_tissue, -array_row, -array_col, -location) %>% 
    tidyr::gather("type", "marker", -pxl_col, -pxl_row)
  return(.data)
}


#'
#'
#'
switch_locations <- function(.data){
  rownames(.data) <- tp[rownames(.data), "location"]
  return(.data)
}


#' Impute, filter and normalize
#' Given a matrix of counts, this will 
#' impute gene counts using knn_smoothing
#' filter based on criteria that gene has
#' 2 reads in atleast 10 cells
#' and finally, quantile normalize the matrix
#' @param counts_mat a dataframe of counts, cells are columns, genes are rows
#' @param qc do some filtering ?
#' @param normalize do one of 'scale, scale-quantile, quantile, sct' normalization, default: None
#' @param only_hvg return only HVG from SCTransform
#' @param ... arguments to pass to k_smooth like k, d and seed
preprocess_expr_mat <- function(counts_mat, impute = T, k = 10, d = 10, seed = 42, 
                                doQC = F, normalize = NULL, only_hvg = T){ 

     if (impute){
         if(!exists('knn_smoothing'))
             stop('make sure you have sourced the imputation file')
         message('smoothing with ', k, ' neighbours and ', d, ' dimensions ...')
         counts_mat <- knn_smoothing(as.matrix(counts_mat), 
                                     k = k, d = d,
                                     seed = seed)
     }
  
  if (doQC) {
    message('performing filtering ...')
    genes <- which(rowSums(counts_mat > 2) > 10) 
    counts_mat <- counts_mat[genes, ]
  }
  
  if (!is.null(normalize)){
    stopifnot(normalize %in% c('scale', 'quantile', 'scale-quantile', 'sct'))
    message('performing ', normalize, ' normalization ...')
    counts_mat <- normalize.decon(counts_mat, 
                                  doQC = F,  
                                  method = normalize, only_hvg = only_hvg)
  }
  
  return(counts_mat)  
  
}


#' normalize a genes by cell matrix
#' optionally perfrom some QC
#' @param data a genes x cell matrix
#' @param doQC should QC be performed: min.cells = 10%, min.features = 500
#' @param method choose a method from 'scale', 'quantile', 'scale-quantile', 'sct
normalize.decon <- function(data, 
                            doQC = FALSE, 
                            method = c('scale', 'quantile', 'scale-quantile', 'sct'), 
                            only_hvg = T){
  
  method = match.arg(method)
  if (method %in% c('scale', 'scale-quantile'))
    scale.factor = median(colSums(data))
  
  if (doQC){
    data <- Seurat::CreateSeuratObject(data, 
                                       min.cells = 0.1*ncol(data), 
                                       min.features = 500, 
                                       assay = 'RNA')
  } else {
    data <- Seurat::CreateSeuratObject(data, 
                                       min.cells = 0, 
                                       min.features = 0, 
                                       assay = 'RNA')  
  }
  
  if (method == 'scale'){
    data <- Seurat::NormalizeData(data, 
                                  normalization.method = 'RC', 
                                  scale.factor = scale.factor)
    data <- as.data.frame(Seurat::GetAssayData(data, 
                                               assay = 'RNA', 
                                               slot = 'data')) 
  } else if (method == 'quantile'){
    data <- as.data.frame(Seurat::GetAssayData(data, 
                                               assay = 'RNA', 
                                               slot = 'counts'))
    data <- normalize.quantiles2(data)
  } else if (method == 'scale-quantile'){
    data <- Seurat::NormalizeData(data, 
                                  normalization.method = 'RC', 
                                  scale.factor = scale.factor)
    data <- as.data.frame(Seurat::GetAssayData(data, 
                                               assay = 'RNA', 
                                               slot = 'data'))
    data <- normalize.quantiles2(data)
  } else if (method == 'sct') {
    message('remember that sct performs some filtering')
    data <- Seurat::SCTransform(data, assay = 'RNA', return.only.var.genes = only_hvg)
    data <- as.data.frame(Seurat::GetAssayData(data, 
                                               assay = 'SCT', 
                                               slot = 'scale.data'))
  }
  
  return(data)
}


#' Normalize quantiles wrapper
#' to handle loss of rownames
#' and colnames
#' @param counts_mat a matrix of couts to be normalized
normalize.quantiles2 <- function(counts_mat){
  RNs <- rownames(counts_mat)
  CNs <- colnames(counts_mat)
  convert = F
  if (!is.matrix(counts_mat)) {
    convert = T
    counts_mat <- as.matrix(counts_mat)
  }
  counts_mat <- preprocessCore::normalize.quantiles(counts_mat)
  rownames(counts_mat) = RNs
  colnames(counts_mat) = CNs
  if (convert) counts_mat <- as.data.frame(counts_mat)
  return(counts_mat)
}


#' make signature matrix from
#' imputed, filtered, normalized 
#' single cell counts
#' @param sc_counts a matrix/dataframe of single cell counts, with cells as columns and genes as rows
#' @param sc_labels a matrix/dataframe of cell type labels, two cols, one for cell and one for cell type
#' @param save_file optionally save the signature matrix
make_signature_matrix <- function(sc_counts, sc_labels, compute_std = F, save_prefix = NULL){
  
  if (is.matrix(sc_counts)) sc_counts <- as.data.frame(sc_counts)
  if (is.matrix(sc_labels)) sc_labels <- as.data.frame(sc_labels)
  colnames(sc_labels) <- c('cell', 'bio_celltype')
  
  sig_mat <- sc_counts %>% 
    t %>% as.data.frame %>%
    tibble::rownames_to_column('cell') %>%
    dplyr::inner_join(sc_labels, by = 'cell') %>% 
    dplyr::select(bio_celltype, everything()) %>%
    dplyr::select(-cell)
  
  sig_mat_means <- sig_mat %>%  
    dplyr::group_by(bio_celltype) %>%
    dplyr::summarise_all(mean) %>%
    tibble::column_to_rownames('bio_celltype') %>%
    t %>% as.data.frame() 
  
  sig_mat_std <- NULL
  if (compute_std){
    sig_mat_std <- sig_mat %>%  
      dplyr::group_by(bio_celltype) %>%
      dplyr::summarise_all(sd) %>%
      tibble::column_to_rownames('bio_celltype') %>%
      t %>% as.data.frame()
  }
  
  if(!is.null(save_prefix)){
    saveRDS(sig_mat_means, file = paste0(save_prefix, '_means.rds'))
    if (compute_std)
      saveRDS(sig_mat_std, file = paste0(save_prefix, '_std.rds'))
  }
  
  invisible(list(sig_means = sig_mat_means, sig_std = sig_mat_std))
}


#' visualization of spots expression
#' across tissue for a single variable
#' @param df spot by value matrix, spots are rows, values are cols
#' @param tp matrix of spot locations, should contain
#'           pxl_col, pxl_row indication column and row location of pixel
#' @param trans any transofrmation that needs to be done to values - log10 etc.
#' @param legend.label label for the value legend
plot.spots0 <- function(df, tp, trans = 'identity', legend.label = 'marker exp.', size_m = 4){
  p <- df %>%
    merge_with_locations(tp) %>%
    ggplot(aes(pxl_row, -pxl_col, color = marker)) +
    geom_point(size = size_m) +
    scale_colour_gradient(low="#f8d39f", high = "#800000", trans = trans) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 14),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.background = element_rect(fill = "black")) +
    labs(colour = legend.label)
  return(p)
}


#' visualization of spots expression
#' across tissue
#' @param df spot by value matrix, spots are rows, values are cols
#' @param tp matrix of spot locations, should contain
#'           pxl_col, pxl_row indication column and row location of pixel
#' @param trans any transofrmation that needs to be done to values - log10 etc.
#' @param legend.label label for the value legend
#' @param joint_scale plot common scale for all spatial distributions
plot.spots <- function(df, tp, trans = 'identity', 
                       legend.label = 'marker exp.', ncols = 3, joint_scale = FALSE, my_theme = NULL, size_m = 1){
  default_theme <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #strip.text = element_text(size = 14),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "black"))
  if (is.null(my_theme)) my_theme = default_theme
  df <- df %>%
    merge_with_locations(tp) 
  if (!joint_scale){
    p <-  df %>%
      group_split(type) %>%
      purrr::map(
        ~ggplot(., aes(pxl_row, -pxl_col, color = marker)) +
          geom_point(size = size_m, shape=16) +
          scale_colour_gradient(low="#f8d39f", high = "#800000", trans = trans) +
          facet_wrap(~ type, labeller = function(x) label_value(x, multi_line = FALSE)) +
          labs(colour = legend.label) +
          my_theme
      ) %>%
      cowplot::plot_grid(plotlist = ., align = 'hv', ncol = ncols)
  } else {
    p <- df %>% ggplot(aes(pxl_row, -pxl_col, color = marker)) +
      geom_point(size = size_m, shape=16) +
      scale_colour_gradient(low="#f8d39f", high = "#800000", trans = trans) +
      facet_wrap(~ type, labeller = function(x) label_value(x, multi_line = FALSE), ncol = ncols) +
      my_theme +
      labs(colour = legend.label)
  }
  return(p)
}


#' plotting correlations of estimate fracs with true values for different celltypes
#' @param results_file an RDS file with the results matrix
#' @param estimates_name ylabel for the plot
#' @param truefractions truth values for the proportions
#' @param col_name column name for error column (estimate - true) default: error
plot.corrs <- function(results_file, estimate_names, truefractions, col_name = 'error'){
  estimates <- readRDS(results_file) %>% as.data.frame() %>%
    tibble::rownames_to_column("spot") %>%
    tidyr::gather("type", "estimates", -spot) %>% 
    inner_join(truefractions, by = c("spot", "type")) %>%
    mutate(!!col_name := (estimates - truefraction)) %>%
    tibble::as_tibble()
  
  p <- ggplot(estimates, aes(x = truefraction, y = estimates)) +
    facet_wrap(~type, scales = "fixed", ncol = 3) +
    geom_point() +
    ggpubr::stat_cor() + geom_abline(slope = 1, intercept = 0, colour = "red") +
    ylab(estimate_names) + xlab("True Proportions")
  print(p)
  invisible(estimates)
}


#' plotting correlations of estimated fracs from different methods
#' @param estimates_file1 an RDS file with the estimates1 matrix
#' @param estimate1_name name for estimates1
#' @param estimates_file2 an RDS file with the estimates2 matrix
#' @param estimate2_name name for estimates2
plot.estimate.corrs <- function(estimates_file1,  estimates_file2, 
                                estimate1_name = 'estimate1', 
                                estimate2_name = 'estimate2'){
  estimates1 <- readRDS(estimates_file1) %>% as.data.frame() %>%
    tibble::rownames_to_column("spot") %>%
    tidyr::gather("type", "estimates1", -spot)
  p <- readRDS(estimates_file2) %>% as.data.frame() %>%
    tibble::rownames_to_column("spot") %>%
    tidyr::gather("type", "estimates2", -spot) %>%
    inner_join(estimates1, by = c("spot", "type")) %>%
    ggplot(aes(x = estimates1, y = estimates2)) +
    facet_wrap(~type, scales = "fixed", ncol = 3) +
    geom_point() +
    ggpubr::stat_cor(method='spearman', cor.coef.name = 'rho') +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    xlab(estimate1_name) + ylab(estimate2_name)
  print(p)
}


#' computes correlation between ST counts and signature counts
#' and plots them
#' @param st_counts gene counts for spots
#' @param sig_mat signatures for different cell types
estimate_plot_corrs <- function(st_counts, sig_mat){
  
  genes <- intersect(rownames(st_counts), rownames(sig_mat))
  df <- apply(st_counts[genes, ], 2, function(x) cor(x, sig_mat[genes, ])) %>%
    t %>% as.data.frame()
  colnames(df) <- colnames(sig_mat)
  
  p <- df %>% tidyr::gather(key = 'type', value = 'corr') %>%
    ggplot(aes(x = corr, fill = type, color = type)) + 
    geom_density(alpha = .3) +
    theme(strip.text = element_text(size = 14),
          title = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16))
  p
}


#' Find average marker gene expression per spot from st_counts
#' @param st_counts matrix of genes x spots counts
#' @param markers output from Seurat::FindAllMarkers()
#' @param topHits number of top markers to consider
#' @param cluster_names optionally supply names for the tSNE/UMAP clusters 
get_marker_matrix <- function(st_counts, markers, topHits = 10, cluster_names = NULL){
  st_counts <- st_counts %>% t %>% as.data.frame()
  # get marker genes
  markers <- markers %>% 
    group_by(cluster) %>% 
    slice_min(p_val_adj, n = topHits) %>% 
    select(gene) %>% as_tibble()
  # make stmarkerexpression matrix
  stMarkerExpression = data.frame(row.names = rownames(st_counts))
  for (celltype in levels(markers$cluster)){
    genes <- markers %>% filter(cluster == celltype) %>% select(gene) %>% unlist
    genes <- intersect(genes, names(st_counts))
    stMarkerExpression[celltype] <- rowMeans(st_counts[genes]) 
  }
  
  return(stMarkerExpression)
}


#' renormalize the proportions matrix
#' using a predefined threshold
#' @param st_props matrix of spot x celltype proportions
#' @param threshold what proportions should be shrunk to zero
renorm_props <- function(st_props, threshold = 0.05){
  st_props[st_props < threshold] <- 0.0
  out <- apply(st_props, 1, function(x) x/sum(x)) 
  out <- t(out)
  if(is.data.frame(st_props)){
    colnames(out) <- colnames(st_props)
    rownames(out) <- rownames(st_props)
  }
}


#' read in MEX format data produced by spaceranger
#' source: 
#' https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/matrices
#' @param matrix_dir directory containing counts, barcode and features
load10xMatrix <- function(matrix_dir) {
  barcode.path <- file.path(matrix_dir, 'barcodes.tsv')
  features.path <- file.path(matrix_dir, 'features.tsv')
  matrix.path <- file.path(matrix_dir, 'matrix.mtx')
  mat <- Matrix::readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE, 
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V2
  return(mat)
}

#' read in the tissue positions file produced by spaceranger
#' ensures consistent columnn nomenclature
#' @param tp_file tissue positions fils from spaceranger
load10xTP <- function(tp_file) {
  columns <- c('barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_col', 'pxl_row')
  tp <- read.csv(tp_file, 
                 col.names = columns, 
                 row.names = 1)
  return(tp)
}


#'
#' @param til_fname
loadTILmap <- function(til_fname){
  read.table(til_fname, col.names = c('x', 'y', 'p', 'blah'))
}


#' Helper function computing overlap between tp square and til square
overlap <- function(tp, til){
  left_x  = max(tp['x'], til['x'])
  right_x = min(tp['x'] + tp['width'], til['x'] + til['width'])
  
  top_y = min(tp['y'] + tp['height'], til['y'] + til['height'])
  bottom_y = max(tp['y'], til['y'])
  
  if (left_x > right_x | bottom_y > top_y)
    return(0)
  
  pxls = (right_x - left_x) * (top_y - bottom_y)
  return(pxls)
}


#' Map TIL values to ST data
#' @param tp st tissue positions file
#' @param til quip prediction file
#' @param spot_size pixel size overlapping a spot (equivalent of square overlapping spot)
#' @param quip_size pixel size of prediction tile (equivalent of 50 microns)
map_til <- function(tp, til, spot_size = 100, quip_size = 50, threshold = 0.094, save_fname = "../data/til"){
  # pxl_row corresponds to the x-coordinate
  # pxl_col corresponds to the y-coordinate
  
  tp <- tp %>% tibble::rownames_to_column('barcode') %>%
    mutate(x = pxl_row - spot_size/2,
           y = pxl_col - spot_size/2,
           height = spot_size,
           width = spot_size) %>%
    tibble::column_to_rownames('barcode') %>%
    filter(in_tissue == 1)
  
  til <- til %>% 
    mutate(height = quip_size,
           width = quip_size,
           p = 0.95*(p > threshold) + (p < threshold)*0.001)
  
  ## set up future
  future::plan(future::multiprocess, workers = 16)

  # should generate a matrix of weights
  # each column is a barcode
  # each row is the overlap with til map box
  weights <- future.apply::future_apply(tp, 1, function(tp_row){
    future.apply::future_apply(til, 1, function(til_row){
      overlap(tp_row, til_row)
    })
  })
  colnames(weights) = rownames(tp)
  saveRDS(weights, file = paste0(save_fname, "_weights.rds"))
  
  out <- apply(weights, 2, function(w){
    if (!all(w == 0))
      weighted.mean(til$p, w)
    else 
      0
  })
  
  out <-  as.data.frame(out)
  saveRDS(out, file = paste0(save_fname, "_prob.rds"))
  out
}

my_sample <- function(.data, num=400){
  if (nrow(.data) < num)
    return(.data)
  dplyr::sample_n(.data, num)
}

