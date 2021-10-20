plot_lbls <- function(lbls, tp){
  lbls %>% merge_with_locations(tp) %>% 
    ggplot(mapping = aes(pxl_row, -pxl_col, color = marker)) +
    geom_point(size = 5) +
    scale_color_brewer(palette="Set2") +
    facet_wrap(~type) +
    theme(
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.background = element_rect(fill = "white")) + 
    labs(colour = 'annotated ct') 
}


#'
#' load the pathology labels for Her2st slides
loadHer2stLbl <- function(lbl_f) {
  read.table(lbl_f, header = T, row.names = 1, sep = "\t", stringsAsFactors = F) %>%
    filter(!is.na(x)) %>% 
    mutate(barcode = paste0(round(x), "x", round(y))) %>%
    tibble::rownames_to_column("temp") %>% 
    select(-temp) %>% 
    tibble::column_to_rownames("barcode") %>% 
    select(label)
}


#' loads the tissue positions for her2st files
#' @param tp_file tissue positions for her2st files
loadHer2stTP <- function(tp_file){
  read.table(tp_file, header = T) %>% 
    mutate(barcode = paste0(x, "x", y)) %>% 
    tibble::column_to_rownames("barcode") %>% 
    dplyr::rename(pxl_row = pixel_x,
           pxl_col = pixel_y,
           in_tissue = selected)
}


#' reads in the proportion estimates rds file and outputs the matrix of proportions
#' in addition, the following columns are output:
#' cancer = sum_i(epithelial_i)
#' tcell.bcell = Tcell + Bcell
#' tcell.bcell.macrophage = Tcell + Bcell + macrophage
#' REMOVED: label : this column contains a label based on proportions, one of c("invasive cancer", "connective tissue", "immune infiltrate")
#' in addition, the column epithelial_4 is renamed to epithelial.4
loadHer2stProp <- function(st_prop){
  if (is.character(st_prop))
    df = readRDS(st_prop)
  else 
    df = st_prop
  #lookup <- c("invasive cancer", "connective tissue", "immune infiltrate")
  df %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("barcode") %>% 
    mutate(cancer = epithelial_1 + epithelial_2 + epithelial_3 + epithelial_4 + epithelial_5,
           tcell.bcell = Tcell + Bcell, 
           tcell.bcell.macrophage = Tcell + Bcell + macrophage) %>%
    #rowwise() %>% mutate(label = lookup[which.max(c(epithelial_4, stroma, tcell.bcell))]) %>% 
    rename(epithelial.4 = epithelial_4) %>%
    tibble::column_to_rownames("barcode")
}


#' 
make_dot_plot <- function(st_props, path_lbl, plt_title){
  p <- merge(st_props, path_lbl, by = 0) %>%
    tibble::column_to_rownames('Row.names') %>% 
    #    select(-base_lbl, -cancer) %>% 
    group_by(path_lbl) %>% 
    summarise_all(median) %>%
    tidyr::gather(key = "cell_type", value = "median_prop", -path_lbl) %>% 
    ggplot(aes(x = path_lbl, y = cell_type)) + 
    geom_point(aes(colour = median_prop, size = median_prop)) + 
    theme_bw(20) +
    rotateTextX() + 
    ggtitle(plt_title)
  print(p)
}


#'
make_dot_plot2 <- function(tbl){
  ggplot(tbl, aes(x = type, y = label)) + 
    geom_point(aes(colour = Freq, size = Freq)) +
    rotateTextX()
  }


#' st_til should have column names "out" for tcell_bcell probabilities
compare_prior_to_baseline <- function(st_props, st_til, path_lbl, plt_title) {
  df <- reduce(list(st_til, st_props, path_lbl), ~merge(.x, .y, by=0) %>% 
                 tibble::column_to_rownames("Row.names")) %>% 
    rename(dl_prob = out, 
           decon_prop = tcell.bcell) %>% 
    select(dl_prob, decon_prop, path_lbl) %>% 
    tibble::rownames_to_column("barcode") %>% 
    tidyr::gather(key = "model", value = "prop", -path_lbl, -barcode)
  
  p1 <- ggplot(df, aes(x=model, y = prop)) +
    geom_boxplot() +
    theme_bw(20) + 
    facet_wrap(~path_lbl, scales = "free_y") + 
    ggtitle(plt_title)
  print(p1)
  
  p2 <- ggplot(df, aes(x=model, y = prop)) +
    geom_line(aes(group = barcode)) +
    theme_bw(20) + 
    facet_wrap(~path_lbl, scales = "free_y") + 
    ggtitle(plt_title)
  print(p2)
}


#' visualize TIL map classifications
visualize_prior_ratio_metrics <- function(st_til, path_lbl, plt_title, filter = F) {
  merge(st_til, path_lbl, by = 0) %>% 
    tibble::column_to_rownames('Row.names') %>% 
    rename(tcell.bcell = out) %>%
    mutate(spot = ifelse(path_lbl == "immune infiltrate", "immune", "rest")) %>% 
    select(-path_lbl) %>% 
    group_by(spot) %>% 
    summarise_all(list(mean = mean, median = median)) %>% 
    tibble::column_to_rownames("spot") %>% 
    t %>% as.data.frame %>% 
    tibble::rownames_to_column("metric") %>% 
    rowwise %>% 
    mutate(ratio = immune/rest) %>% 
    tidyr::gather(key = type, value, -metric) %>% 
    ggplot(aes(x = type, y = value, fill = metric)) + 
    geom_bar(stat = 'identity', position = position_dodge2()) +
    coord_flip() + theme_bw(20) +
    geom_text(aes(label = round(value, 4)), size = 10, position=position_dodge(width=0.9)) +
    geom_abline(intercept = 1, slope = 0, 
                color = 'red', size = 1.5, 
                linetype = 2) +
    ggtitle(plt_title)
}


visualize_prior_ratio_metrics_cancer <- function(st_til, path_lbl, plt_title, filter = F) {
  merge(st_til, path_lbl, by = 0) %>% 
    tibble::column_to_rownames('Row.names') %>% 
    rename(tcell.bcell = out) %>%
    mutate(spot = ifelse(path_lbl %in% c("invasive cancer", "cancer in situ") , "cancer", "rest")) %>% 
    select(-path_lbl) %>% 
    group_by(spot) %>% 
    summarise_all(list(mean = mean, median = median)) %>% 
    tibble::column_to_rownames("spot") %>% 
    t %>% as.data.frame %>% 
    tibble::rownames_to_column("metric") %>% 
    rowwise %>% 
    mutate(ratio = cancer/rest) %>% 
    tidyr::gather(key = type, value, -metric) %>% 
    ggplot(aes(x = type, y = value, fill = metric)) + 
    geom_bar(stat = 'identity', position = position_dodge2()) +
    coord_flip() + theme_bw(20) +
    geom_text(aes(label = round(value, 4)), size = 10, position=position_dodge(width=0.9)) +
    geom_abline(intercept = 1, slope = 0, 
                color = 'red', size = 1.5, 
                linetype = 2) +
    ggtitle(plt_title)
}


#' 
get_median_props <- function(st_props, path_lbl){
  merge(st_props, path_lbl, by = 0) %>%
    tibble::column_to_rownames('Row.names') %>% 
    select(-prior_lbl, -cancer) %>% 
    group_by(path_lbl) %>% 
    summarise_all(median) %>% 
    filter(path_lbl == "immune infiltrate") %>%
    select(Bcell, macrophage, Tcell, tcell_bcell, tcell_bcell_macrophage)
}


#' summarise proportions by performance metrics
get_perf_metrics <- function(st_props, path_lbl) {
  merge(st_props, path_lbl, by = 0) %>%
    tibble::column_to_rownames('Row.names') %>%
    group_by(path_lbl) %>%
    summarise_all(list(median = median, mean = mean)) %>% 
    filter(path_lbl == "immune infiltrate") %>% 
    select(Bcell_median, macrophage_median, Tcell_median, tcell.bcell_median, tcell.bcell.macrophage_median,
           Bcell_mean, macrophage_mean, Tcell_mean, tcell.bcell_mean, tcell.bcell.macrophage_mean, 
           epithelial.4_mean, epithelial.4_median, cancer_mean, cancer_median) 
}


get_perf_metrics_cancer <- function(st_props, path_lbl) {
  merge(st_props, path_lbl, by = 0) %>%
    mutate(spot = ifelse(path_lbl %in% c("invasive cancer", "cancer in situ") , "cancer", "rest")) %>%
    select(-path_lbl) %>%
    tibble::column_to_rownames('Row.names') %>%
    group_by(spot) %>%
    summarise_all(list(median = median, mean = mean)) %>% 
    filter(spot == "cancer") %>% 
    select(Bcell_median, macrophage_median, Tcell_median, tcell.bcell_median, tcell.bcell.macrophage_median,
           Bcell_mean, macrophage_mean, Tcell_mean, tcell.bcell_mean, tcell.bcell.macrophage_mean, 
           epithelial.4_mean, epithelial.4_median, cancer_mean, cancer_median) 
}


get_ratio_metric <- function (st_props, path_lbl, filter = F, spot_type = c("immune", "cancer")) {
  st_props <- merge(st_props, path_lbl, by = 0) %>%
    tibble::column_to_rownames('Row.names')
  spot_type = match.arg(spot_type)
  if (spot_type == "immune") 
    cell_types = c("immune infiltrate")
  else if (spot_type == "cancer")
    cell_types = c("invasive cancer", "cancer in situ")

  if (filter)
    st_props <- st_props %>% filter(path_lbl != "undetermined")
  
  st_props %>%
    mutate(spot = ifelse(path_lbl %in% cell_types, spot_type, "rest")) %>% 
    select(-path_lbl) %>%
    group_by(spot) %>% 
    summarise_all(list(median = median, mean = mean))  %>% 
    select(spot, 
           Bcell_median, macrophage_median, Tcell_median, tcell.bcell_median, tcell.bcell.macrophage_median,
           Bcell_mean, macrophage_mean, Tcell_mean, tcell.bcell_mean, tcell.bcell.macrophage_mean, 
           epithelial.4_mean, epithelial.4_median, cancer_mean, cancer_median) %>% 
    tibble::column_to_rownames("spot") %>%
    t %>% data.frame() %>% 
    tibble::rownames_to_column("metric") %>%
    tidyr::separate(metric, c("celltype", "summary"), sep = "_") %>%
    rowwise %>%
    mutate(ratio = (!!sym(spot_type))/rest)  
}


get_ratio_metric2 <- function (st_props, path_lbl, filter = F) {
  st_props <- merge(st_props, path_lbl, by = 0) %>%
    tibble::column_to_rownames('Row.names') 
  
  if (filter)
    st_props <- st_props %>% filter(path_lbl != "undetermined")
  
  st_props %>%
    mutate(spot = ifelse(path_lbl == "immune infiltrate", "immune", "rest")) %>% 
    select(-path_lbl) %>%
    group_by(spot) %>% 
    summarise_all(list(median = median, mean = mean))  %>% 
    select(spot, 
           Bcell_median, macrophage_median, Tcell_median, tcell.bcell_median, tcell.bcell.macrophage_median,
           Bcell_mean, macrophage_mean, Tcell_mean, tcell.bcell_mean, tcell.bcell.macrophage_mean, 
           epithelial.4_mean, epithelial.4_median, cancer_mean, cancer_median) %>% 
    tibble::column_to_rownames("spot") %>%
    t %>% data.frame() %>% 
    tibble::rownames_to_column("metric") %>%
    rowwise %>%
    mutate(ratio = immune/rest ) %>%
    tibble::column_to_rownames("metric") %>% 
    select(ratio)  %>% 
    t %>% as.data.frame
}


get_metrics = function(st_props, path_lbl, metric_function, name = "base") {
  df = metric_function(st_props, path_lbl) %>% 
    t %>% as.data.frame %>% 
    rename(!!rlang::sym(paste0(name, "_ratio")) := ratio) %>% 
    merge(get_perf_metrics(st_props, path_lbl) %>% 
            t %>% as.data.frame %>% 
            rename(!!rlang::sym(paste0(name, "_summary")) := V1), by = 0) %>% 
    rename(metric = Row.names) %>% 
    tidyr::separate(metric, c("celltype", "metric"), sep = "_")
  return(df)
}


visualize_metrics <- function(st_props, path_lbl, title, metric_func) {
  st_metrics <- metric_func(st_props, path_lbl) %>% 
    t %>% as.data.frame %>%
    tibble::rownames_to_column("metric") %>% 
    tidyr::separate(metric, c("cell_type", "summary"), extra = "merge", sep = "_")
  colnames(st_metrics) <- c("cell_type", "summary", "value")
  p <- st_metrics %>% 
    ggplot(aes(x = cell_type, y = value, fill = summary)) + 
    geom_bar(stat = 'identity', position = position_dodge2()) +
    coord_flip() + theme_bw(20) + 
    geom_abline(intercept = 1, slope = 0, 
                color = 'red', size = 1.5, 
                linetype = 2) +
    ggtitle(title) +
    geom_text(aes(label = round(value, 2)), size = 15, position=position_dodge(width=0.9))
  # st_metrics <- st_metrics %>% 
  #   rename(baseline = value) %>% 
  #   tidyr::separate(metric, c("cell_type", "summary"), extra = "merge", sep = "_")
  print(p)
  invisible(st_metrics %>% rename(baseline = value))
}


plot_metric_by_lambda <- function (df, metric_func, baseline_df, plt_title) {
  props <- map(df, ~metric_func(.x, path_lbl))
  sapply(props, rbind) %>% 
    as.data.frame %>%
    tibble::rownames_to_column("cell_type") %>%
    tidyr::separate(cell_type, c("cell_type", "summary"), extra = "merge", sep = "_") %>%
    left_join(baseline_df, by = c("cell_type", "summary")) %>%
    tidyr::gather(key = "model", value = "prop", -cell_type, -summary, -baseline) %>% 
    mutate(prior_lambda = as.numeric(stringr::str_replace(model, "model_lambda_", "")),
           prop = as.numeric(prop)) %>%
    ggplot(aes(x = prior_lambda, y = prop, color = summary)) + 
    geom_line() + 
    facet_wrap(~cell_type, scales = "free_y") + 
    theme_bw(20) + ggtitle(plt_title) +
    geom_hline(aes(yintercept = baseline, colour = summary), linetype = 2, size = 2)
}

read_dl_pred_ptch <- function(pred_f){
  read.table(pred_f, col.names = c("x", "y", "p", "blah"))
}

plot_dl_ptch_on_img <- function(img, dl_pred_ptch, ...){
  info <- image_info(img)
  ggplot(dl_pred_ptch, aes(x = x, y = -y, color = p, alpha = ifelse(p > 0.1, 2, p)) ) +
    annotation_raster(img, 0, info$width, 0, -info$height) + 
    geom_point(...) +
    scale_colour_gradient(low = "#f8d39f", high = "#800000") + 
    theme_void()
}


plot_dl_pred_on_img <- function(img, dl_pred, tp, size_m = 5) {
  info <- image_info(img)
  dl_pred_tp <- merge(dl_pred, tp, by = 0) %>% tibble::column_to_rownames("Row.names") %>% rename(til_prob = out)
  ggplot(dl_pred_tp, aes(x = pxl_row, y = -pxl_col, color = til_prob, alpha = til_prob)) + 
    annotation_raster(img, 0, info$width, 0, -info$height) + 
    geom_point(size = size_m) +
    scale_colour_gradient(low = "#f8d39f", high = "#800000") + 
    theme_void()
}


plot_decon_pred_on_img <- function(img, decon_pred, tp, size_m = 5) {
  info <- image_info(img)
  decon_pred_tp <- merge(decon_pred, tp, by = 0) %>% tibble::column_to_rownames("Row.names") 
  ggplot(decon_pred_tp, aes(x = pxl_row, y = -pxl_col, color = tcell.bcell, alpha = tcell.bcell)) + 
    annotation_raster(img, 0, info$width, 0, -info$height) + 
    geom_point(size = size_m) +
    scale_colour_gradient(low = "#f8d39f", high = "#800000") + 
    theme_void()
}


plot_lbls_on_img <- function(img, path_lbl, tp, size_m = 5, ...){
  info <- image_info(img)
  path_lbl_tp <- merge(path_lbl, tp, by = 0) %>% tibble::column_to_rownames("Row.names")
  ggplot(path_lbl_tp, aes(x = pxl_row, y = -pxl_col, color = label)) + 
    annotation_raster(img, 0, info$width, 0, -info$height) + 
    geom_point(size = size_m, ...) +
    theme_void() + 
    # previous colours: low: #f8d39f  
    scale_color_manual(values = c("Non-Immune" = "#fdae6b", "Immune" = "#800000"))
}



#' slides: entries are filenames 
#' slides: names are slide names
combine_slides = function(slides){
  modify_counts = function(x, y){
    cnts = read.table(file.path("~/projects/her2st/data/ST-cnts", x))
    rownames(cnts) = paste0(y, "_", rownames(cnts))
    cnts %>% t %>% as.data.frame()
  }
  slide_cnts = imap(slides, ~modify_counts(.x, .y))
  reduce(slide_cnts, ~merge(.x, .y, by = 0) %>% tibble::column_to_rownames("Row.names"))
}


#' mega_slide (barcodes x cell_types)
uncombine_slides = function(prefixes, mega_slide){
  barcodes = rownames(mega_slide)
  results = list()
  for (p in prefixes){
    row_names = stringr::str_starts(barcodes, p)
    sub_s = mega_slide[row_names, ]
    rownames(sub_s) = gsub(rownames(sub_s), pattern = paste0(p, "_"), replacement = "")
    results[[p]] = sub_s
  }
  results
}


her2st2tenx <- function(tp){
  tp %>% 
    tibble::rownames_to_column("barcode") %>%
    dplyr::rename(array_row = new_x) %>%
    dplyr::rename(array_col = new_y) %>%
    select(barcode, in_tissue, array_row, array_col, pxl_col, pxl_row)
}


move_images <- function(from_loc, to_loc, barcodes){
  if (!dir.exists(to_loc))
    dir.create(to_loc)
  for (b in barcodes){
    fname = paste0(b, ".png")
    from = file.path(from_loc, fname)
    to = file.path(to_loc, fname)
    file.copy(from, to)
  }
}


create_dir <- function(path_lbls_f, st_props_f, from_loc, to_loc, random = c("boxplot", "ringo"), save_fname = "path_images.rds"){
  random = match.arg(random)
  
  ## move annotated images
  path_lbl = loadHer2stLbl(path_lbls_f)
  barcodes_anno = path_lbl %>% 
    dplyr::filter(label == "immune infiltrate") %>%
    tibble::rownames_to_column("barcode") %>% 
    select(barcode) %>% 
    pull
  move_images(from_loc, to_loc, barcodes_anno)
  
  ## move outliers
  st_props = loadHer2stProp(st_props_f)
  if (random == "boxplot") {
    out <- boxplot(st_props$tcell.bcell, main = "distribution of TIL proportions", col = "blue", cex =1.5, pch = 16)
    out <- out$stats[5]
    barcodes <- st_props %>% 
      filter(tcell.bcell > out) %>% 
      tibble::rownames_to_column("barcode") %>% 
      select(barcode) %>%
      pull
    barcodes_out = setdiff(barcodes, barcodes_anno)
    move_images(from_loc, to_loc, barcodes_out)
  } else {
    out <- Ringo::upperBoundNull(st_props$tcell.bcell, limits = c(0,1))
    hist(st_props$tcell.bcell, 30, col = "blue", main = "Histogram of TIL proportions", xlab = "")
    abline(v = out, col = "red", lwd = 2, lty = 2)
    barcodes <- st_props %>% 
      filter(tcell.bcell > out) %>% 
      tibble::rownames_to_column("barcode") %>% 
      select(barcode) %>%
      pull
    barcodes_out = setdiff(barcodes, barcodes_anno)
    move_images(from_loc, to_loc, barcodes_out)
  }
  
  ## move random
  remaining_set = setdiff(rownames(path_lbl), c(barcodes_anno, barcodes_out))
  n_out = length(barcodes_out)
  if (! length(remaining_set) < n_out){
    barcodes_random = sample(remaining_set, n_out)
    move_images(from_loc, to_loc, barcodes_random)
  } else{
    cat("not enough remaining barcodes ... \n")
    barcodes_random = NULL
  }

  n_bar_anno = length(barcodes_anno)
  n_bar_out = length(barcodes_out)
  n_bar_random = length(barcodes_random)
  
  df = data.frame(barcode = c(barcodes_anno, barcodes_out, barcodes_random),
                  type = c(rep("annotated", n_bar_anno), rep("outlier", n_bar_out), rep("random", n_bar_random)))
  saveRDS(df, file = file.path(to_loc, save_fname))

  cat("moved", n_bar_anno, "annotated images,", 
      "moved", n_bar_out, "outlier images,", 
      "moved", n_bar_random, "random images\n")
}


anonymize <- function(patch_dir, anon_fname = "anonymized_fnames.rds"){
  l = list.files(patch_dir, pattern = "*.png")
  s = sample(seq(length(l)), replace = F)
  from_files = file.path(patch_dir, l)
  to_files = file.path(patch_dir, paste0("image", "_", s, ".png"))
  purrr::walk2(from_files, to_files, ~file.rename( .x, .y))
  df = data.frame(orig_fname = l,
                  anon_name = paste0("image", "_", s, ".png"))
  saveRDS(df, file = file.path(patch_dir, anon_fname))
}


unanonymize <- function(scores_dir, anon_fname){
  df = readRDS(anon_fname) %>% tibble::column_to_rownames("anon_name")
  dirs = c("low", "medium", "high")
  out = list()
  f = function(x) file.path(scores_dir, d, x)
  for (d in dirs){
    ddf  = df[list.files(file.path(scores_dir, d), "*.png"), , drop = F] %>% tibble::rownames_to_column("anon_name")
    with(ddf, purrr::walk2(orig_fname, anon_name, ~file.rename(f(.y), f(.x))))
    out[[d]] = ddf$orig_fname
  }
  out = purrr::imap(out, ~data.frame(barcode = unlist(strsplit(.x, ".png")), label = .y)) %>% dplyr::bind_rows()
  out
}


compare_low <- function(scores_dir, anon_fname, algo_fname){
  df = readRDS(anon_fname) %>% tibble::column_to_rownames("anon_name")
  none_dir = file.path(scores_dir, "low", "none")
  none_imgs = list.files(none_dir, "*.png")
  df = df[none_imgs, , drop = F] %>%
    mutate(label = "none") %>%
    rename(barcode = orig_fname) %>% rowwise() %>%
    mutate(barcode = strsplit(barcode, ".png")[[1]])
  p = readRDS(algo_fname)
  merge(df, p, by = "barcode") %>% 
    select(label, type) %>% 
    table
}
 

#' get p-values using a permuation scheme
get_p_values <- function(gistEsts, exprEsts, pathologistImmuneSpots, n_perms = 100000, summary_func = c(mean, median)) {
  ratioGist0 <- summary_func(gistEsts[pathologistImmuneSpots == 1]) / summary_func(gistEsts[pathologistImmuneSpots == 0])
  ratioExpr0 <- summary_func(exprEsts[pathologistImmuneSpots == 1]) / summary_func(exprEsts[pathologistImmuneSpots == 0])
  observedDiff <- ratioGist0 - ratioExpr0
  
  permDiff <- numeric()
  for(i in seq(n_perms)){
    pathologistImmuneSpotsPerm <- sample(pathologistImmuneSpots)
    ratioGist <- summary_func(gistEsts[pathologistImmuneSpotsPerm == 1]) / summary_func(gistEsts[pathologistImmuneSpotsPerm == 0])
    ratioExpr <- summary_func(exprEsts[pathologistImmuneSpotsPerm == 1]) / summary_func(exprEsts[pathologistImmuneSpotsPerm == 0])
    permDiff[i] <- ratioGist - ratioExpr
  }
  
  pval = (sum(permDiff > abs(observedDiff)) + sum(permDiff < -abs(observedDiff))) / n_perms
  
  invisible(list(permDiff = permDiff,
                 ratioGist = ratioGist0, 
                 ratioExpr = ratioExpr0, 
                 observedDiff = observedDiff, 
                 pval = pval))  
}
