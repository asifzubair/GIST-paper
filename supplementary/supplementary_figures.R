############################################################
#
# This script generates the supplementary figures, 
# including all the panel figures. Files are numbered
# according to figure's occurence in the supplementary text
#
############################################################


## libraries
library(purrr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(ggExtra)
library(ggcorrplot)
library(imager)
library(magick)

## helper functions
source("../utils/st_helpers.R")
source("../utils/tmb2_helpers.R")
source("../utils/her2st_helpers.R")

if (!dir.exists("./export"))
  dir.create("./export")

#### Suppl. Fig. 1

## Data
results <- readRDS("./data/results-splatter-plus.29Apr2020.rds")
cols <- paste0("Group", seq(6))
rows <- paste0("Bulk", seq(100))
results <- lapply(results[1:8], function(x){colnames(x) <- cols; rownames(x) <- rows; x})

gold_standard <- results$truth
truefractions <- gold_standard %>% as.data.frame() %>% tibble::rownames_to_column("sample") %>% 
  tidyr::gather("type", "truefraction", -sample)

types = c(Group1 = "Cell type 1",
          Group2 = "Cell type 2",
          Group3 = "Cell type 3",
          Group4 = "Cell type 4",
          Group5 = "Cell type 5",
          Group6 = "Cell type 6")

linreg <- results$lm %>% as.data.frame() %>%
  tibble::rownames_to_column("sample") %>% 
  tidyr::gather("type", "estimates", -sample) %>% 
  inner_join(truefractions, by=c("sample", "type")) %>%
  mutate(lm = (estimates - truefraction))

cbs <- results$cbs %>% as.data.frame() %>%
  tibble::rownames_to_column("sample") %>% 
  tidyr::gather("type", "estimates", -sample) %>% 
  inner_join(truefractions, by=c("sample", "type")) %>%
  mutate(Cibersort = (estimates - truefraction))

drnas <- results$drnas %>% as.data.frame() %>%
  tibble::rownames_to_column("sample") %>% 
  tidyr::gather("type", "estimates", -sample) %>% 
  inner_join(truefractions, by=c("sample", "type")) %>%
  mutate(DeconRNAseq = (estimates - truefraction))

estimates <- readr::read_tsv('data/W.2020-06-20154131.618764.tsv') %>% 
  dplyr::rename(sample = X1)

stereo <- estimates %>% 
  tidyr::gather("type", "estimates", -sample) %>% 
  inner_join(truefractions, by=c("sample", "type")) %>%
  mutate(Stereoscope = (estimates - truefraction))

out <- results$bayconMean %>% as.data.frame() %>%
  tibble::rownames_to_column("sample") %>% 
  tidyr::gather("type", "estimates", -sample) %>% 
  inner_join(truefractions, by=c("sample", "type")) %>%
  mutate(`GIST-base` = (estimates - truefraction))

## Figure

## (a)
p <- ggplot(linreg, aes(x=truefraction, y=estimates)) +
  facet_wrap(~ type, scales = "fixed", ncol = 3, labeller = as_labeller(types)) +
  geom_point(size = 0.02) +
  xlab("True proportions") + 
  ylab("Linear regression estimates") + 
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  geom_abline(slope = 1, intercept = 0, colour = "red", size = 0.2) +
  stat_cor(size = 1.25, method = "spearman") + 
  theme(strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text = element_text(family="Arial", face="plain", size=8), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.ticks = element_line(colour = "black", size = 0.1),
        axis.line = element_line(colour="black", size = 0.3),
        axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
        axis.title = element_text(family="Arial", face="plain", size=8)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  ggExtra::rotateTextX()
ggsave(filename="export/fig01a_lin_reg.svg", 
       plot=p, height=3, width=3, units=c("in"), dpi=300)

## (b)
p <- ggplot(cbs, aes(x=truefraction, y=estimates)) +
  facet_wrap(~ type, scales = "fixed", ncol = 3, labeller = as_labeller(types)) +
  geom_point(size = 0.02) +
  xlab("True proportions") + 
  ylab("CIBERSORT estimates") + 
  ylim(0,1) + 
  xlim(0,1) + 
  geom_abline(slope = 1, intercept = 0, colour = "red", size = 0.2) +
  stat_cor(size = 1, method = "spearman") + 
  theme(strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text = element_text(family="Arial", face="plain", size=8), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.ticks = element_line(colour = "black", size = 0.1),
        axis.line = element_line(colour="black", size = 0.3),
        axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
        axis.title = element_text(family="Arial", face="plain", size=8)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  ggExtra::rotateTextX()
ggsave(filename="export/fig01b_cbs.svg", 
       plot=p, height=3, width=3, units=c("in"), dpi=300)


## (c)
p <- ggplot(drnas, aes(x=truefraction, y=estimates)) +
  facet_wrap(~ type, scales = "fixed", ncol = 3, labeller = as_labeller(types)) +
  xlab("True proportions") + 
  ylab("DeconRNAseq estimates") + 
  ylim(0,1) + 
  xlim(0,1) + 
  geom_abline(slope = 1, intercept = 0, colour = "red", size = 0.2) +
  geom_point(size = 0.02) +
  stat_cor(size = 1, method = "spearman") + 
  theme(strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text = element_text(family="Arial", face="plain", size=8), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.ticks = element_line(colour = "black", size = 0.1),
        axis.line = element_line(colour="black", size = 0.3),
        axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
        axis.title = element_text(family="Arial", face="plain", size=8)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  ggExtra::rotateTextX()
ggsave(filename="export/fig01c_drnas.svg", 
       plot=p, height=3, width=3, units=c("in"), dpi=300)  

## (d)
p <- ggplot(stereo, aes(x=truefraction, y=estimates)) +
  facet_wrap(~ type, scales = "fixed", ncol = 3, labeller = as_labeller(types)) +
  xlab("True proportions") + 
  ylab("Stereoscope estimates") + 
  ylim(0,1) + 
  xlim(0,1) + 
  geom_abline(slope = 1, intercept = 0, colour = "red", size = 0.2) +
  geom_point(size = 0.02) +
  stat_cor(size = 1, method = "spearman") + 
  theme(strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text = element_text(family="Arial", face="plain", size=8), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.ticks = element_line(colour = "black", size = 0.1),
        axis.line = element_line(colour="black", size = 0.3),
        axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
        axis.title = element_text(family="Arial", face="plain", size=8)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  ggExtra::rotateTextX()
ggsave(filename="export/fig01d_stereo.svg", 
       plot=p, height=3, width=3, units=c("in"), dpi=300)

## (e)
p <- ggplot(out, aes(x=truefraction, y=estimates)) +
  facet_wrap(~ type, scales = "fixed", ncol = 3, labeller = as_labeller(types)) +
  xlab("True proportions") + 
  ylab("GIST-base estimates") + 
  ylim(0,1) + 
  xlim(0,1) + 
  geom_abline(slope = 1, intercept = 0, colour = "red", size = 0.2) +
  geom_point(size = 0.02) +
  stat_cor(size = 1, method = "spearman") + 
  theme(strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text = element_text(family="Arial", face="plain", size=8), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.ticks = element_line(colour = "black", size = 0.1),
        axis.line = element_line(colour="black", size = 0.3),
        axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
        axis.title = element_text(family="Arial", face="plain", size=8)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  ggExtra::rotateTextX()
ggsave(filename="export/fig01e_gist.svg", 
       plot=p, height=3, width=3, units=c("in"), dpi=300)

#### Suppl. Fig. 2

## Data
load("./data/idb_mixing_study_mapped_results.RData")
cols <- c("value", paste0("value", seq(99)))

truefractions <- gold_standard %>% tibble::rownames_to_column("type") %>% 
  tidyr::gather("sample", "truefraction", -type)
colnames(linreg) <- cols
colnames(drnas) <- cols
colnames(out) <- cols

linreg <- linreg %>% tibble::rownames_to_column("type") %>% 
  tidyr::gather("sample", "estimates", -type) %>% 
  inner_join(truefractions, by=c("sample", "type")) %>%
  mutate(lm = (estimates - truefraction))

cbs <- cbs %>% tibble::rownames_to_column("type") %>% 
  tidyr::gather("sample", "estimates", -type) %>% 
  inner_join(truefractions, by=c("sample", "type")) %>%
  mutate(Cibersort = (estimates - truefraction))

drnas <- drnas %>% tibble::rownames_to_column("type") %>% 
  tidyr::gather("sample", "estimates", -type) %>% 
  inner_join(truefractions, by=c("sample", "type")) %>%
  mutate(DeconRNAseq = (estimates - truefraction))

out <- out %>% tibble::rownames_to_column("type") %>% 
  tidyr::gather("sample", "estimates", -type) %>% 
  inner_join(truefractions, by=c("sample", "type")) %>%
  mutate(`GIST-base` = (estimates - truefraction))

## Figure

## (a)
p <- ggplot(linreg, aes(x=truefraction, y=estimates)) +
  facet_wrap(~ type, scales = "fixed", ncol = 3) +
  xlab("True proportions") + 
  ylab("Linear regression estimates") + 
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  geom_abline(slope = 1, intercept = 0, colour = "red", size = 0.2) +
  geom_point(size = 0.02) +
  stat_cor(size = 1.25, method = "spearman") + 
  theme(strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text = element_text(family="Arial", face="plain", size=6), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.ticks = element_line(colour = "black", size = 0.1),
        axis.line = element_line(colour="black", size = 0.3),
        axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
        axis.title = element_text(family="Arial", face="plain", size=8)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  ggExtra::rotateTextX()
ggsave(filename="export/fig02a_lin_reg.svg", 
       plot=p, height=3.5, width=4, units=c("in"), dpi=300)

## (b)
p <- ggplot(cbs, aes(x=truefraction, y=estimates)) +
  facet_wrap(~ type, scales = "fixed", ncol = 3) +
  xlab("True proportions") + 
  ylab("CIBERSORT estimates") + 
  ylim(0,1) + 
  xlim(0,1) + 
  geom_abline(slope = 1, intercept = 0, colour = "red", size = 0.2) +
  geom_point(size = 0.02) +
  stat_cor(size = 1.25, method = "spearman") + 
  theme(strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text = element_text(family="Arial", face="plain", size=6), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.ticks = element_line(colour = "black", size = 0.1),
        axis.line = element_line(colour="black", size = 0.3),
        axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
        axis.title = element_text(family="Arial", face="plain", size=8)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  ggExtra::rotateTextX()
ggsave(filename="export/fig02b_cbs.svg", 
       plot=p, height=3.5, width=4, units=c("in"), dpi=300)


## (c)
p <- ggplot(drnas, aes(x=truefraction, y=estimates)) +
  facet_wrap(~ type, scales = "fixed", ncol = 3) +
  xlab("True proportions") + 
  ylab("DeconRNAseq estimates") + 
  ylim(0,1) + 
  xlim(0,1) + 
  geom_abline(slope = 1, intercept = 0, colour = "red", size = 0.2) +
  geom_point(size = 0.02) +
  stat_cor(size = 1.25, method = "spearman") + 
  theme(strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text = element_text(family="Arial", face="plain", size=6), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.ticks = element_line(colour = "black", size = 0.1),
        axis.line = element_line(colour="black", size = 0.3),
        axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
        axis.title = element_text(family="Arial", face="plain", size=8)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  ggExtra::rotateTextX()
ggsave(filename="export/fig02c_drnas.svg", 
       plot=p, height=3.5, width=4, units=c("in"), dpi=300)  

## (d)
p <- ggplot(out, aes(x=truefraction, y=estimates)) +
  facet_wrap(~ type, scales = "fixed", ncol = 3) +
  xlab("True proportions") + 
  ylab("GIST-base estimates") + 
  ylim(0,1) + 
  xlim(0,1) + 
  geom_abline(slope = 1, intercept = 0, colour = "red", size = 0.2) +
  geom_point(size = 0.02) +
  stat_cor(size = 1.25, method = "spearman") + 
  theme(strip.background = element_rect(fill = "white", colour = "white"), 
        strip.text = element_text(family="Arial", face="plain", size=6), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.ticks = element_line(colour = "black", size = 0.1),
        axis.line = element_line(colour="black", size = 0.3),
        axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
        axis.title = element_text(family="Arial", face="plain", size=8)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  ggExtra::rotateTextX()
ggsave(filename="export/fig02d_gist.svg", 
       plot=p, height=3.5, width=4, units=c("in"), dpi=300)


#### Suppl. Fig. 3

## Export theme
my_theme = theme(legend.key.width = unit(0.06, "in"),
                 legend.key.height = unit(0.09, "in"),
                 legend.key = element_rect(colour = "white", fill = "white"),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 4),
                 legend.box.background = element_rect(colour = "white"),
                 legend.box.margin = margin(0,0,0,0), 
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_rect(fill="white", colour="white"),
                 strip.text = element_text(size = 6),
                 axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.ticks = element_blank(),
                 panel.background = element_rect(fill = "black"),
                 plot.margin = unit(c(0, 0, 0, 0), "cm"))

## Data
st_props_base = readRDS("data/results_baycon_base_decon_10xMouseBrain2_sig_mat_stereo_sc_k5d10_sct")
tp <- load10xTP("data/tissue_positions_list.csv")

## Figure
p <- plot.spots(st_props_base, tp, legend.label = 'Prop.', ncols = 6, my_theme = my_theme, size_m = 0.1)
ggsave(filename = "export/fig03_gist.svg", p, height = 11, width = 8, units = c("in"), dpi = 300)


#### Suppl. Fig. 4

## Data
estimatedProportions = as.matrix(readRDS("data/estimatedProportions_RCTD.rds"))  

## Figure
p <- plot.spots(estimatedProportions, tp, legend.label = 'Prop.', ncols = 6, my_theme = my_theme, size_m = 0.1)
ggsave(filename = "export/fig04_rctd.svg", p, height = 11, width = 8, units = c("in"), dpi = 300)


#### Suppl. Fig. 5

## Data
spotlight_ls <- readRDS("data/spotlight_ls.rds")
decon_mtrx <- spotlight_ls[[2]]
rownames(decon_mtrx) = readRDS("data/st_counts_barcodes.rds")
estimatedProportions = decon_mtrx
colnames(estimatedProportions) <- gsub("\\.", "_", colnames(estimatedProportions))

## Figure
p <- plot.spots(estimatedProportions[, -57], tp, legend.label = 'Prop.', ncols = 6, my_theme = my_theme, size_m = 0.1)
ggsave(filename = "export/fig05_spotlight.svg", p, height = 11, width = 8, units = c("in"), dpi = 300)


#### Suppl. Fig. 6

## Data
cell2loc_file = "data/W_cell_density.csv"
estimatedProportions = read.csv(cell2loc_file, header = 1) %>% 
  mutate(spot_id = gsub("tenx_mouse_brain2_", "", spot_id)) %>%
  tibble::column_to_rownames("spot_id") %>%
  rename_all(function(x) gsub("mean_spot_factors", "", x))

## Figure
p <- plot.spots(estimatedProportions, tp, legend.label = 'Prop.', ncols = 6, my_theme = my_theme, size_m = 0.1)
ggsave(filename = "export/fig06_cell2location.svg", p, height = 11, width = 8, units = c("in"), dpi = 300)


#### Suppl. Fig. 7

## Data
estimatedProportions <- read.table('data/W.2020-12-01143927.020572.tsv')

## Figure
p <- plot.spots(estimatedProportions, tp, legend.label = 'Prop.', ncols = 6, my_theme = my_theme, size_m = 0.1)
ggsave(filename = "export/fig07_stereo.svg", p, height = 11, width = 8, units = c("in"), dpi = 300)


#### Suppl. Fig. 8 to Suppl. Fig. 12

## Data
fignums = c("08", "09", "10", "11", "12")
slides = c("A1", "D1", "E1", "F1", "H1")
img_pths = c("data/9769_C1_HE_small_edited_transparent_outline_black.jpg",
             "data/10426_HE_BT_C1_edited_transparent_outline_black.jpg",
             "data/HE_BT23567_D2_edited_transparent_outline_black.jpg",
             "data/HE_BT23810_D2_edited_transparent_outline_black.jpg",
             "data/HE_BT24044_D2_edited_transparent_outline_black.jpg") 
imgs = map(img_pths, image_read)

path_gt_pths = paste0("data/", slides, "_labeled_coordinates.tsv")
path_gts = map(path_gt_pths, loadHer2stLbl)

tp_img_pths = paste0("data/", slides , "_selection.tsv") 
tp_imgs = map(tp_img_pths, loadHer2stTP)

pred_fs = paste0("data/prediction-", slides)
dl_pred_ptchs <- map(pred_fs, read_dl_pred_ptch)

dl_pred_img_pths = paste0("data/", tolower(slides), "_til_mapped_decon_k10d10_sct.rds")
dl_pred_imgs = map(dl_pred_img_pths, readRDS)

decon_img_pths = paste0("data/results_mega_slide_", tolower(slides), "_k10d10_sct.rds")
decon_imgs = map(decon_img_pths, loadHer2stProp)

decon_gist_pths = paste0("data/results_baycon_base_decon_stdata_mega_slide_", tolower(slides),
                         "_k10d10_noFilter_sct_beta_prior_mapped_probabilities_lambda_50_tnbc_sigmat_sct.rds")
decon_gists = map(decon_gist_pths, loadHer2stProp)

pval_pths <- paste0("data/", tolower(slides), "_pval_median_summary.rds")
pvals <- map(pval_pths, readRDS)

## Figure

## Suppl. Fig. 8b to Suppl. Fig. 12b
for (i in seq_along(slides)){
  fignum = fignums[i]
  s = slides[i]
  img = imgs[[i]]
  tp_img = tp_imgs[[i]]
  path_gt = path_gts[[i]] 
  path_gt <- path_gt %>% 
    tibble::rownames_to_column("spot") %>%
    mutate(label = ifelse(label == "immune infiltrate", "Immune", "Non-Immune")) %>% 
    tibble::column_to_rownames("spot") 
  
  p <- plot_lbls_on_img(img, path_gt, tp_img, 0.6, alpha = 1) +
    theme(
      legend.box.margin = margin(5.25,14,5.25,14),
      legend.key.width = unit(0.08, "in"),
      legend.key.height = unit(0.06, "in"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(family="Arial", face="plain", size=8),
      legend.title = element_blank()
    ) +
    guides(size = FALSE, alpha = FALSE) + 
    labs(color = "Annotation")
  
  ggsave(filename=paste0("export/fig", fignum, "b_", s, "_path_anno.png"), 
         plot=p, height=2, width=2.3, units=c("in"), dpi=600)
}


## Suppl. Fig. 8c to Suppl. Fig. 12c
for (i in seq_along(slides)){
  fignum = fignums[i]
  s = slides[i]
  img = imgs[[i]]
  tp_img = tp_imgs[[i]]
  pxl_row_max =  max(tp_img$pxl_row)
  pxl_row_min = min(tp_img$pxl_row)
  pxl_col_max = max(tp_img$pxl_col)
  pxl_col_min = min(tp_img$pxl_col)
  
  dl_pred_ptch <- dl_pred_ptchs[[i]]
  dl_pred_ptch <- dl_pred_ptch %>% 
    filter(x >= max(0, (pxl_row_min - 200))) %>% 
    filter(x <= (pxl_row_max + 100)) %>% 
    filter(y >= max(0, (pxl_col_min - 200))) %>% 
    filter(y <= (pxl_col_max + 100))
  p <- plot_dl_ptch_on_img(img, dl_pred_ptch, shape = 15, size = 0.25) +
    theme(
      legend.box.margin = margin(5.25,14,5.25,14),
      legend.key.width = unit(0.25, "in"),
      legend.key.height = unit(0.06, "in"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(family="Arial", face="plain", size=8),
      legend.title = element_blank(),
      title = element_text(family="Arial", face="plain", size=8),
      plot.title = element_text(hjust = 0.5),
    ) +
    guides(size = FALSE, alpha = FALSE)
  ggsave(filename=paste0("export/fig", fignum, "c_", s, "_ptch_til", ".png"), 
         plot=p, height=2.15, width=2.3, units=c("in"), dpi=600)
}


## Suppl. Fig. 8d to Suppl. Fig. 12d
for (i in seq_along(slides)){
  fignum = fignums[i]
  s = slides[i]
  img = imgs[[i]]
  dl_pred_img = dl_pred_imgs[[i]]
  tp_img = tp_imgs[[i]]
  p <- plot_dl_pred_on_img(img, dl_pred_img, tp_img, 0.6) +
    theme(
      legend.box.margin = margin(5.25,14,5.25,14),
      legend.key.width = unit(0.25, "in"),
      legend.key.height = unit(0.06, "in"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(family="Arial", face="plain", size=8),
      legend.title = element_blank(),
      title = element_text(family="Arial", face="plain", size=8),
      plot.title = element_text(hjust = 0.5),
    ) +
    guides(size = FALSE, alpha = FALSE)
  ggsave(filename = paste0("export/fig", fignum, "d_", s, "_til", ".png"), plot=p, height=2.15, width=2.3, units=c("in"), dpi=600)
}


## Suppl. Fig. 8e to Suppl. Fig. 12e
upper_lims = vector("double", length(slides))
for (i in seq_along(slides)){
  fignum = fignums[i]
  s = slides[i]
  img = imgs[[i]]
  decon_img = decon_imgs[[i]]
  tp_img = tp_imgs[[i]]
  upper_lims[i] = decon_imgs[[i]] %>% select(tcell.bcell) %>% max
  upper_lim = upper_lims[i]
  p <- plot_decon_pred_on_img(img, decon_img, tp_img, 0.6) +
    theme(
      legend.box.margin = margin(5.25,14,5.25,14),
      legend.key.width = unit(0.25, "in"),
      legend.key.height = unit(0.06, "in"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(family="Arial", face="plain", size=8),
      legend.title = element_blank(),
      title = element_text(family="Arial", face="plain", size=8),
      plot.title = element_text(hjust = 0.5),
    ) +
    guides(size = FALSE, alpha = FALSE) + 
    scale_colour_gradient(low="#f8d39f", high="#800000", limits=c(0, upper_lim))
  ggsave(filename = paste0("export/fig", fignum, "e_", s, "_base_decon", ".png"), 
         plot=p, height=2.15, width=2.3, units=c("in"), dpi=600)
}


## Suppl. Fig. 8f to Suppl. Fig. 12f
for (i in seq_along(slides)){
  fignum = fignums[i]
  s = slides[i]
  img = imgs[[i]]
  decon_gist = decon_gists[[i]]
  tp_img = tp_imgs[[i]]
  upper_lim = upper_lims[i]
  p <- plot_decon_pred_on_img(img, decon_gist, tp_img, 0.6) +
    theme(
      legend.box.margin = margin(5.25,14,5.25,14),
      legend.key.width = unit(0.25, "in"),
      legend.key.height = unit(0.06, "in"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(family="Arial", face="plain", size=8),
      legend.title = element_blank(),
      title = element_text(family="Arial", face="plain", size=8),
      plot.title = element_text(hjust = 0.5),
    ) +
    guides(size = FALSE, alpha = FALSE) + 
    scale_colour_gradient(low="#f8d39f", high="#800000", limits=c(0, upper_lim))
  ggsave(filename = paste0("export/fig", fignum, "f_", s, "_gist", ".png"), 
         plot=p, height=2.15, width=2.3, units=c("in"), dpi=600)  
}


## Suppl. Fig. 8g to Suppl. Fig. 12g
for (i in seq_along(slides)){
  fignum = fignums[i]
  s = slides[i]
  o <- pvals[[i]]
  p <- ggplot(mapping = aes(o$permDiff)) +
    geom_histogram(bins = 100, color = "darkblue") + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.background = element_rect(fill = "white", colour = "white"),
      axis.ticks = element_line(colour = "black", size = 0.1),
      axis.line = element_line(colour="black", size = 0.3),
      axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
      axis.title = element_text(family="Arial", face="plain", size=8)) +
    xlab("Permuted test statistic") +
    ylab("Count") +
    geom_vline(xintercept = o$observedDiff, color = "red", size = 0.2)
  ggsave(filename = paste0("export/fig", fignum, "g_", s, "_pval", ".svg"), 
         plot=p, height=2, width=2.3, units=c("in"), dpi=300)
}


## Suppl. Fig. 8h to Suppl. Fig. 12h
for (i in seq_along(slides)){
  fignum = fignums[i]
  s = slides[i]
  dl_pred_img = dl_pred_imgs[[i]]
  decon_base_img = decon_imgs[[i]]
  p <- merge(dl_pred_img, decon_base_img %>% select(tcell.bcell), by = 0) %>% 
    tibble::column_to_rownames("Row.names") %>% 
    ggplot(aes(x = tcell.bcell, y = out)) +
    geom_point(size = 0.1) +
    ggpubr::stat_cor(method='spearman', cor.coef.name = 'rho', size = 3, 
                     label.x.npc = 'left', label.y.npc = 'top') +
    geom_smooth(method = "lm", level=0, color = "red", size = 0.2) + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = "white"),
      axis.ticks = element_line(colour = "black", size = 0.1),
      axis.line = element_line(colour="black", size = 0.3),
      axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
      axis.title = element_text(family="Arial", face="plain", size=8)) +
    xlab("Gene expression-derived\n immune cell proportions") +
    ylab("Deep learning-derived\n immune cell proportions")
  ggsave(filename = paste0("export/fig", fignum, "h_", s, "_corr", ".svg"), 
         plot=p, height=2, width=2.3, units=c("in"), dpi=300)
}


#### Suppl. Fig. 13

## Data
decon_prior_b1 = loadHer2stProp(
  "data/results_baycon_base_decon_stdata_mega_slide_b1_k10d10_noFilter_sct_beta_prior_mapped_probabilities_lambda_50_tnbc_sigmat_sct.rds")
o_b1 = readRDS("data/path_annotation_B1.rds")
a_b1 = readRDS("data/path_images_B1.rds")
df_b1 = purrr::reduce(list(decon_prior_b1 %>% tibble::rownames_to_column("barcode"), o_b1, a_b1), function(x, y) merge(x, y, by = "barcode"))

img_c1 = image_read("data/5714_HE_BT_C1_edited_transparent_outline_black.jpg")
decon_base_c1 = loadHer2stProp("data/results_mega_slide_c1_k10d10_sct.rds")
tp_c1 = loadHer2stTP("data/C1_selection.tsv")

dl_pred_c1 = readRDS("data/c1_til_mapped_decon_k10d10_sct.rds")

decon_prior_c1 = loadHer2stProp(
  "data/results_baycon_base_decon_stdata_mega_slide_c1_k10d10_noFilter_sct_beta_prior_mapped_probabilities_lambda_50_tnbc_sigmat_sct.rds")

table_c1 = readRDS("data/table_C1.rds") 
colnames(table_c1) <- c("GIST", "Random")
rownames(table_c1) <- c("Low", "Middle", "High")
table_c1 = table_c1 %>% as.data.frame()

o_c1 = readRDS("data/path_annotation_C1.rds")
a_c1 = readRDS("data/path_images_C1.rds")
df_c1 = purrr::reduce(list(decon_prior_c1 %>% tibble::rownames_to_column("barcode"), o_c1, a_c1), function(x, y) merge(x, y, by = "barcode"))

## Figure

## (a)
p <- df_b1 %>% 
  ggplot(aes(y = tcell.bcell, 
             x = forcats::fct_recode(forcats::fct_relevel(label, c("low", "medium", "high")), Low="low", Middle="medium", High="high"), 
             fill = type)) +
  geom_boxplot(position = position_dodge(1), lwd = 0.1, outlier.size = 0.01) +
  xlab("path_score") +
  theme(legend.key.size = unit(0.1, 'in'),
        legend.key = element_rect(colour = "white", fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(family="Arial", face="plain", size=8),
        legend.box.background = element_rect(colour = "white"),
        legend.box.margin = margin(0,0,0,0),
        legend.position = "bottom",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.ticks = element_line(colour = "black", size = 0.1),
        axis.line = element_line(colour="black", size = 0.3),
        axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
        axis.title = element_text(family="Arial", face="plain", size=8)) + 
  ylab("GIST score") +
  xlab("2nd pathologist's score") +
  scale_fill_manual(values = c("#fff7bc", "#99d8c9"), labels = c("GIST", "Random")) +
  ggExtra::rotateTextX()
ggsave(filename="export/fig13a_gist_score_annotation_b1.svg", 
       plot=p, height=2.3, width=2.25, units=c("in"), dpi=300)

## (b)
p <- plot_decon_pred_on_img(img_c1, decon_base_c1, tp_c1, 0.6) +
  theme(
    legend.box.margin = margin(5.25,14,5.25,14),
    legend.key.width = unit(0.25, "in"),
    legend.key.height = unit(0.06, "in"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(family="Arial", face="plain", size=8),
    legend.title = element_blank(),
    title = element_text(family="Arial", face="plain", size=8),
    plot.title = element_text(hjust = 0.5),
  ) +
  guides(size = FALSE, alpha = FALSE) 
ggsave(filename="export/fig13b_gist_base_on_c1.png", 
       plot=p, height=2.15, width=2.3, units=c("in"), dpi=600)

## (c)
p <- plot_dl_pred_on_img(img_c1, dl_pred_c1, tp_c1, 0.6) +
  theme(
    legend.box.margin = margin(5.25,14,5.25,14),
    legend.key.width = unit(0.25, "in"),
    legend.key.height = unit(0.06, "in"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(family="Arial", face="plain", size=8),
    legend.title = element_blank(),
    title = element_text(family="Arial", face="plain", size=8),
    plot.title = element_text(hjust = 0.5),
  ) +
  guides(size = FALSE, alpha = FALSE) 
ggsave(filename="export/fig13c_til_on_c1.png", 
       plot=p, height=2.15, width=2.3, units=c("in"), dpi=600)

## (d)
p <- merge(dl_pred_c1, decon_base_c1 %>% select(tcell.bcell), by = 0) %>% 
  tibble::column_to_rownames("Row.names") %>% 
  ggplot(aes(x = tcell.bcell, y = out)) +
  geom_point(size = 0.1) +
  ggpubr::stat_cor(method='spearman', cor.coef.name = 'rho', size = 3, 
                   label.x.npc = 'left', label.y.npc = 'top') +
  geom_smooth(method = "lm", level=0, color = "red", size = 0.2) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "white", colour = "white"),
    axis.ticks = element_line(colour = "black", size = 0.1),
    axis.line = element_line(colour="black", size = 0.3),
    axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
    axis.title = element_text(family="Arial", face="plain", size=8)) +
  xlab("Gene expression-derived\n immune cell proportions")+ 
  ylab("Deep learning-derived\n immune cell proportions")
ggsave(filename="export/fig13d_corr_til_decon_on_c1.svg", 
       plot=p, height=2, width=2.3, units=c("in"), dpi=300)

## (e)
p <- plot_decon_pred_on_img(img_c1, decon_prior_c1, tp_c1, 0.6) +
  theme(
    legend.box.margin = margin(5.25,14,5.25,14),
    legend.key.width = unit(0.25, "in"),
    legend.key.height = unit(0.06, "in"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(family="Arial", face="plain", size=8),
    legend.title = element_blank(),
    title = element_text(family="Arial", face="plain", size=8),
    plot.title = element_text(hjust = 0.5),
  ) +
  guides(size = FALSE, alpha = FALSE) 
ggsave(filename="export/fig13e_gist_on_c1.png", 
       plot=p, height=2.15, width=2.3, units=c("in"), dpi=600)

## (f)
p <- make_dot_plot2(table_c1) + 
  guides(size=F) +
  theme(legend.key.size = unit(0.1, 'in'),
        legend.key = element_rect(colour = "white", fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(family="Arial", face="plain", size=8),
        legend.box.background = element_rect(colour = "white"),
        legend.box.margin = margin(0,0,0,0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.ticks = element_line(colour = "black", size = 0.1),
        axis.line = element_line(colour="black", size = 0.3),
        axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
        axis.title.y = element_text(family="Arial", face="plain", size=8),
        axis.title.x = element_blank())+ 
  ylab("2nd pathologist's\n immune cell score") +
  scale_colour_gradient(low="#a8ddb5", high = "#43a2ca")
ggsave(filename="export/fig13f_dot_pot_c1.svg", 
       plot=p, height=2, width=2.25, units=c("in"), dpi=300)

## (g)
p <- df_c1 %>% 
  ggplot(aes(y = tcell.bcell, 
             x = forcats::fct_recode(forcats::fct_relevel(label, c("low", "medium", "high")), Low="low", Middle="medium", High="high"), 
             fill = type)) +
  geom_boxplot(position = position_dodge(1), lwd = 0.1, outlier.size = 0.01) +
  xlab("path_score") +
  theme(legend.key.size = unit(0.1, 'in'),
        legend.key = element_rect(colour = "white", fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(family="Arial", face="plain", size=8),
        legend.box.background = element_rect(colour = "white"),
        legend.box.margin = margin(0,0,0,0),
        legend.position = "bottom",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.ticks = element_line(colour = "black", size = 0.1),
        axis.line = element_line(colour="black", size = 0.3),
        axis.text = element_text(family="Arial", face="plain", size=8, colour="black"),
        axis.title = element_text(family="Arial", face="plain", size=8)) + 
  ylab("GIST score") +
  xlab("2nd pathologist's score") +
  scale_fill_manual(values = c("#fff7bc", "#99d8c9"), labels = c("GIST", "Random")) +
  ggExtra::rotateTextX()
ggsave(filename="export/fig13g_gist_score_annotation_c1.svg", 
       plot=p, height=2.3, width=2.5, units=c("in"), dpi=300)
