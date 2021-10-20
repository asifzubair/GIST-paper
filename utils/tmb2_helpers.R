#neuron_index = 16:45
#glial_index = c(1:5, 7, 46:50)
glial_prefix = c("Astrocytes", "Ependymal", "Oligos")


get_props2 <- function(df){
  df %>% 
    as.data.frame %>%
    tibble::rownames_to_column('spot') %>%
    mutate(neuron = rowSums(select(., starts_with("Neurons"))),
           glial = rowSums(select(., starts_with(glial_prefix)))) %>% 
    select(neuron, glial, spot)  %>% 
    tibble::column_to_rownames('spot')
}


get_props <- function(df){
  df %>% 
    as.data.frame %>%
    tibble::rownames_to_column('spot') %>%
    mutate(neuron = rowSums(select(., starts_with("Neurons"))),
           glial = rowSums(select(., starts_with(glial_prefix)))) %>%
    select(neuron, glial, spot) %>%
    tidyr::gather("type", "estimates", -spot)  
}


iwalk_plot <- function(st_props, label = "ihc", label.x.npc = "center", label.y.npc = "bottom") {
  iwalk(st_props, ~print(ggplot(.x, aes(ihc, estimates)) + 
                           facet_wrap(~type, scales = "fixed", ncol = 3) +
                           geom_point() +
                           ggpubr::stat_cor(method='spearman', cor.coef.name = 'rho', size = 8, 
                                            label.x.npc = label.x.npc, label.y.npc = label.y.npc) +
                           geom_abline(slope = 1, intercept = 0, colour = "red") + 
                           ggtitle(.y) + 
                           theme_bw(24) + 
                           xlab(label)
                       )
      )
}


qq_ihc <- function(st_props, st_ihc, ct, ...){
  x = st_props %>% filter(type == ct) %>% pull(estimates)
  y = st_ihc %>% filter(type == ct) %>% pull(ihc)
  qqplot(x, y, xlab = 'Deconvolution proportions', ylab = 'IHC', ...)
  abline(0, 1, col = 'red')
}


read.VisiumSpatialRNA2 <- function (datadir) {
  coords <- readr::read_csv(file = paste(datadir, "spatial/tissue_positions_list.csv", 
                                         sep = "/"), 
                            col_names = c("barcodes", "in_tissue", "x", "y", "pxl_col_in_fullres", "pxl_row_in_fullres"))
  coords = tibble::column_to_rownames(coords, var = "barcodes") %>% select(x, y)
  counts <- Seurat::Read10X_h5(paste0(datadir, "/filtered_feature_bc_matrix.h5"))
  puck = SpatialRNA(coords, counts)
  restrict_puck(puck, colnames(puck@counts))
}

create.RCTD2 <- function(spatialRNA, reference, max_cores = 8, test_mode = FALSE, gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.75, UMI_min = 100, UMI_max = 200000,
                        class_df = NULL, CELL_MIN_INSTANCE = 25, cell_type_names = NULL) {
  
  config <- list(gene_cutoff = gene_cutoff, fc_cutoff = fc_cutoff, gene_cutoff_reg = gene_cutoff_reg, fc_cutoff_reg = fc_cutoff_reg, UMI_min = UMI_min, max_cores = max_cores,
                 N_epoch = 8, N_X = 50000, K_val = 100, N_fit = 1000, N_epoch_bulk = 30, MIN_CHANGE_BULK = 0.0001, MIN_CHANGE_REG = 0.001, UMI_max = UMI_max, MIN_OBS = 3)
  if(test_mode)
    config <- list(gene_cutoff = .00125, fc_cutoff = 0.5, gene_curoff_reg = 0.002, fc_cutoff_reg = 0.75, UMI_min = 1000,
                   N_epoch = 1, N_X = 50000, K_val = 100, N_fit = 50, N_epoch_bulk = 4, MIN_CHANGE_BULK = 1, MIN_CHANGE_REG = 0.001, UMI_max = 200000, MIN_OBS = 3, max_cores = 1)
  if(is.null(cell_type_names))
    cell_type_names <- levels(reference@cell_types)
  cell_type_info <- list(info = RCTD:::process_cell_type_info(reference, cell_type_names, CELL_MIN = CELL_MIN_INSTANCE), renorm = NULL)
  
  puck = restrict_counts(spatialRNA, rownames(spatialRNA@counts), UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
  print('create.RCTD: getting regression differentially expressed genes: ')
  #puckMeans <- rowMeans(sweep(puck@counts, 2 , puck@nUMI, '/'))
  gene_list_reg = get_de_genes(cell_type_info$info, puck, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg, MIN_OBS = config$MIN_OBS)
  if(length(gene_list_reg) == 0)
    stop("create.RCTD: Error: 0 regression differentially expressed genes found")
  print('create.RCTD: getting platform effect normalization differentially expressed genes: ')
  gene_list_bulk = get_de_genes(cell_type_info$info, puck, fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff, MIN_OBS = config$MIN_OBS)
  if(length(gene_list_bulk) == 0)
    stop("create.RCTD: Error: 0 bulk differentially expressed genes found")
  puck = restrict_counts(puck, gene_list_bulk, UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
  puck = restrict_puck(puck, colnames(puck@counts))
  if(is.null(class_df))
    class_df <- data.frame(cell_type_info$info[[2]], row.names = cell_type_info$info[[2]]); colnames(class_df)[1] = "class"
  internal_vars <- list(gene_list_reg = gene_list_reg, gene_list_bulk = gene_list_bulk, proportions = NULL, class_df = class_df)
  new("RCTD", spatialRNA = puck, reference = reference, config = config, cell_type_info = cell_type_info, internal_vars = internal_vars)
}
