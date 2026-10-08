library(clusterProfiler)
library(tidyverse)
library(cowplot)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
source("02_st-xenium/util.R")
source("utils/plot_gsea_annot.R")
source("utils/run_gsea_cp.R")

## -------Enrichment analysis-------

top_genes <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/xenium_cNMF_k4_topgenes.csv", row.names = 1)
gep_scores <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/xenium_cNMF_k4_gep_scores.csv", row.names = 1)


### -------- GEP4 GO -------------

res_gep4 <- run_go_gsea(setNames(gep_scores$X4, rownames(gep_scores)), species = "human", p_cutoff=0.05)
res_gep4 <- clusterProfiler::simplify(res_gep4)
View(res_gep4@result)

p1 <- plot_gsea_annot(res_gep4, "GO:0045814", color = "#7b0f0f", label_size=4, label_alpha=0.5)
p2 <- plot_gsea_annot(res_gep4, "GO:0006282", color = "#7b0f0f", label_size=4, label_alpha=0.5)
p3 <- plot_gsea_annot(res_gep4, "GO:0006397", color = "#7b0f0f", label_size=4, label_alpha=0.5)

ggsave(plot = p1, path="02_st-xenium/figures", filename= "GO_GSEA_0045814.pdf", width=7, height=5, dpi = 700, units = "cm")
ggsave(plot = p2, path="02_st-xenium/figures", filename= "GO_GSEA_0006282.pdf", width=7, height=5, dpi = 700, units = "cm")
ggsave(plot = p3, path="02_st-xenium/figures", filename= "GO_GSEA_0006397.pdf", width=7, height=5, dpi = 700, units = "cm")


### -------- HALLMARK -------------

res_gep1 <- run_hallmark_gsea(setNames(gep_scores$X1, rownames(gep_scores)), species = "human", p_cutoff=1)
res_gep2 <- run_hallmark_gsea(setNames(gep_scores$X2, rownames(gep_scores)), species = "human", p_cutoff=1)
res_gep3 <- run_hallmark_gsea(setNames(gep_scores$X3, rownames(gep_scores)), species = "human", p_cutoff=1)
res_gep4 <- run_hallmark_gsea(setNames(gep_scores$X4, rownames(gep_scores)), species = "human", p_cutoff=1)

df_all <- bind_rows(
  res_gep1@result %>%
    mutate(GEP = "GEP 1"),
  
  res_gep2@result %>%
    mutate(GEP = "GEP 2"),
  
  res_gep3@result %>%
    mutate(GEP = "GEP 3"),
  
  res_gep4@result %>%
    mutate(GEP = "GEP 4")
)

sig_ids <- df_all %>% filter(p.adjust < 0.05) %>% pull(ID) %>% unique()
nonsig_ids <- setdiff(unique(df_all$ID), sig_ids)
df_all_sig <- df_all %>% filter(!ID %in% nonsig_ids)
df_all_sig$ID <- gsub("HALLMARK_", "", df_all_sig$ID)
NES_M <- df_all_sig %>% dplyr::select(ID, NES, GEP) %>% pivot_wider(names_from = "GEP", values_from = "NES") %>% column_to_rownames("ID")
P_M <- df_all_sig %>% dplyr::select(ID, p.adjust, GEP) %>% pivot_wider(names_from = "GEP", values_from = "p.adjust") %>% column_to_rownames("ID")

head(NES_M)
head(P_M)

#### ----- Dot plot -------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtree)
  library(aplot)
  library(patchwork)
  library(cowplot)
})

# --- Inputs: NES_M and P_M are matrices with identical dimnames ---
stopifnot(all(dim(NES_M) == dim(P_M)))
stopifnot(identical(rownames(NES_M), rownames(P_M)))
stopifnot(identical(colnames(NES_M), colnames(P_M)))

# 1) Cluster pathways by NES pattern across GEPs (correlation distance)
X <- as.matrix(NES_M)

# robust: handle any NAs
cor_mtx <- cor(t(X), use = "pairwise.complete.obs", method = "pearson")
dist_mtx <- as.dist(1 - cor_mtx)        # distance in [0,2]
hc <- hclust(dist_mtx, method = "average")

tree_plot <- ggtree(hc)

# tip order to force dotplot y-axis order
tip_order <- tree_plot$data %>%
  dplyr::filter(isTip) %>%
  dplyr::arrange(y) %>%
  dplyr::pull(label)

# 2) Build long df for dotplot: Description (gene set), GEP (x), NES (color), -log10(p) (size)
df_long <- as.data.frame(NES_M) %>%
  tibble::rownames_to_column("Description") %>%
  tidyr::pivot_longer(-Description, names_to = "GEP", values_to = "NES") %>%
  left_join(
    as.data.frame(P_M) %>%
      tibble::rownames_to_column("Description") %>%
      tidyr::pivot_longer(-Description, names_to = "GEP", values_to = "p_adj"),
    by = c("Description", "GEP")
  ) %>%
  mutate(
    # hide non-significant points
    NES   = ifelse(p_adj < 0.05, NES, NA_real_),
    log10p = ifelse(p_adj < 0.05, -log10(pmax(p_adj, .Machine$double.xmin)), NA_real_),
    Description = factor(Description, levels = rev(tip_order))
  )

# 3) Dotplot styling (same as your template)
sc <- scale_colour_gradient2(low = "blue", mid = "white", high = "red")

dotplot <- ggplot(df_long, aes(x = GEP, y = Description, color = NES, size = log10p)) +
  geom_point() +
  geom_point(shape = 21, colour = "black", stroke = 0.5) +
  cowplot::theme_cowplot() +
  scale_size(range = c(0, 3.5)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  ylab("") + xlab(NULL) +
  scale_y_discrete(position = "right") +
  sc

# 4) Align y-lims and combine tree + dotplot
tree_plot <- tree_plot + aplot::ylim2(dotplot)

p <- tree_plot + dotplot + plot_layout(ncol = 2, widths = c(0.7, 4))
ggsave(plot = p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures/", filename= "GEP_enrichment.pdf", width=13, height=14, dpi = 700, units = "cm")


# GEP 2: ECM.
# GEP 3: Prolif/DDR
# GEP 4: 

# GEP 1: Stress response
#Proteotoxic stress / Heat-shock response:
#  [HSPA6, HSPD1, HSPA9, HSPA4, HSPA8, HSPB8, DNAJA1, DNAJB2, DNAJB4, TCP1, CCT6A, CCT7, CCT5, PTGES3, LONP1, SERPINH1]
#Integrated Stress Response (ISR) / ER stress / UPR (ATF4–CHOP axis):
#  [ATF4, ATF3, DDIT3 (CHOP), XBP1, GADD45A, TRIB3, SESN2, EIF4EBP1]
#Autophagy / Protein degradation:
#  [SQSTM1, MAP1LC3B, UBE2D3, CYLD, RAB7A, RAB5A, PPID, PPIF]

# GEP 2: MES
#ECM remodeling / matrix crosslinking:
#  [LOXL2, LOXL1, LOX, MMP14, TNC, COL7A1, COL5A1, FMOD, P4HA1, SERPINH1, TIMP2, TIMP4]
#TGF-β / mesenchymal transition signaling:
#  [TGFB3, SNAI2, THBS1, ITGA5, ITGA3, CD44, PDPN, TGM2, FSTL1]
#Wnt / Notch / developmental signaling:
# [WNT5A, NOTCH3, DKK3, LGR4]
#Cell–cell / cell–matrix adhesion:
#[LAMB3, ITGA5, ITGA3, CD44, CD82, SDC2, TSPAN9]
#Lipid metabolism / cholesterol biosynthesis:
# [FADS2, SQLE, FASN, EBP, FDFT1, DHCR24, ACLY]

# GEP 3: Inflam.
#Pro-inflammatory cytokines / acute inflammatory response:
#  [IL1B, IL6, IL11, IL24, CSF3, CCL20, BDKRB1, IL1RN]
#	NF-κB / AP-1 immediate early response:
#  [FOSL1, NR4A1, ERRFI1, DUSP6, PDE4B, PTGS2]
#	Growth factor signaling (EGFR / FGFR / HGF axis):
#  [EGFR, EREG, HBEGF, FGF2, FGF7, HGF, SPRY2]
#	Angiogenesis / vascular remodeling:
#  [VEGFA, VEGFC, PTX3, ESM1, EDIL3]

# GEP4: Regulatory
#RNA splicing / mRNA processing:
#  [SF3B1, SRSF1, PRPF31, CPSF6, DDX39B, DDX51, MOV10, PABPC1L, NXF1, SUGP2, LUC7L, ZCCHC8]
#	Chromatin modifiers / epigenetic regulation:
#  [EZH2, DOT1L, PHF8, TRRAP, KAT2A, ASXL1]
#	Transcriptional regulation / nuclear signaling:
#  [STAT2, HMGA2, PVT1, MEG3, CTBP1, IFRD1]

# res_gep1 <- run_c8_gsea(setNames(gep_scores$X1, rownames(gep_scores)), species = "human", p_cutoff=1)
# res_gep2 <- run_c8_gsea(setNames(gep_scores$X2, rownames(gep_scores)), species = "human", p_cutoff=1)
# res_gep3 <- run_c8_gsea(setNames(gep_scores$X3, rownames(gep_scores)), species = "human", p_cutoff=1)
# res_gep4 <- run_c8_gsea(setNames(gep_scores$X4, rownames(gep_scores)), species = "human", p_cutoff=1)


## --------- Patient state genesets ---------

# === Read MP sets ===
sig <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/data/Nomura_2025_signatures.csv")
sig <- sig[,grepl("^MP", colnames(sig))]


### ----- GEP 1 -------

out <- run_gsea_cp(sig, gep_scores, stat_col="X1")
as.data.frame(out$gsea_cp)

p1 <- plot_gsea_annot(out$gsea_cp, "MP_10_Stress1", color = "#7b0f0f", label_size=4, label_alpha=0.5, show_x_axis_elements=FALSE)
p2 <- plot_gsea_annot(out$gsea_cp, "MP_15_Stress2", color = "#7b0f0f", label_size=4, label_alpha=0.5, show_x_axis_elements=FALSE)
p3 <- plot_gsea_annot(out$gsea_cp, "MP_6_MES", color = "#0b3c6f", label_size=4, label_alpha=0.5, show_x_axis_elements=FALSE)

ggsave(plot = p1, path="02_st-xenium/figures", filename= "GEP1_gsea_plot_MP_10_Stress1.pdf", width=6*.9, height=5*.9, dpi = 700, units = "cm")
ggsave(plot = p2, path="02_st-xenium/figures", filename= "GEP1_gsea_plot_MP_15_Stress2.pdf", width=6*.9, height=5*.9, dpi = 700, units = "cm")
ggsave(plot = p3, path="02_st-xenium/figures", filename= "GEP1_gsea_plot_MP_6_MES.pdf", width=6*.9, height=5*.9, dpi = 700, units = "cm")

### ----- GEP 2 -------

out <- run_gsea_cp(sig, gep_scores, stat_col="X2")
as.data.frame(out$gsea_cp)

p1 <- plot_gsea_annot(out$gsea_cp, "MP_6_MES", color = "#7b0f0f", label_size=4, label_alpha=0.5, show_x_axis_elements=FALSE)
p2 <- plot_gsea_annot(out$gsea_cp, "MP_10_Stress1", color = "#0b3c6f", label_size=4, label_alpha=0.5, show_x_axis_elements=FALSE)
p3 <- plot_gsea_annot(out$gsea_cp, "MP_4_AC", color = "#7b0f0f", label_size=4, label_alpha=0.5, show_x_axis_elements=FALSE)

ggsave(plot = p1, path="02_st-xenium/figures", filename= "GEP2_gsea_plot_MP_6_MES.pdf", width=6*.9, height=5*.9, dpi = 700, units = "cm")
ggsave(plot = p2, path="02_st-xenium/figures", filename= "GEP2_gsea_plot_MP_10_Stress1.pdf", width=6*.9, height=5*.9, dpi = 700, units = "cm")
ggsave(plot = p3, path="02_st-xenium/figures", filename= "GEP2_gsea_plot_MP_4_AC.pdf", width=6*.9, height=5*.9, dpi = 700, units = "cm")


### ----- GEP 3 -------

out <- run_gsea_cp(sig, gep_scores, stat_col="X3")
as.data.frame(out$gsea_cp)

p1 <- plot_gsea_annot(out$gsea_cp, "MP_10_Stress1", color = "#0b3c6f", label_size=4, label_alpha=0.5, show_x_axis_elements=FALSE)
p2 <- plot_gsea_annot(out$gsea_cp, "MP_14_NRGN", color = "#7b0f0f", label_size=4, label_alpha=0.5, show_x_axis_elements=FALSE)
p3 <- plot_gsea_annot(out$gsea_cp, "MP_2_OPC", color = "#7b0f0f", label_size=4, label_alpha=0.5, show_x_axis_elements=FALSE)

ggsave(plot = p1, path="02_st-xenium/figures", filename= "GEP3_gsea_plot_MP_6_MES.pdf", width=6*.9, height=5*.9, dpi = 700, units = "cm")
ggsave(plot = p2, path="02_st-xenium/figures", filename= "GEP3_gsea_plot_MP_10_Stress1.pdf", width=6*.9, height=5*.9, dpi = 700, units = "cm")
ggsave(plot = p3, path="02_st-xenium/figures", filename= "GEP3_gsea_plot_MP_4_AC.pdf", width=6*.9, height=5*.9, dpi = 700, units = "cm")


### ----- GEP 4 -------

out <- run_gsea_cp(sig, gep_scores, stat_col="X4")
as.data.frame(out$gsea_cp)

p1 <- plot_gsea_annot(out$gsea_cp, "MP_15_Stress2", color = "#0b3c6f", label_size=4, label_alpha=0.5, show_x_axis_elements=FALSE)
p2 <- plot_gsea_annot(out$gsea_cp, "MP_14_NRGN", color = "#0b3c6f", label_size=4, label_alpha=0.5, show_x_axis_elements=FALSE)
p3 <- plot_gsea_annot(out$gsea_cp, "MP_6_MES", color = "#0b3c6f", label_size=4, label_alpha=0.5, show_x_axis_elements=FALSE)

ggsave(plot = p1, path="02_st-xenium/figures", filename= "GEP4_gsea_plot_MP_15_Stress2.pdf", width=6*.9, height=5*.9, dpi = 700, units = "cm")
ggsave(plot = p2, path="02_st-xenium/figures", filename= "GEP4_gsea_plot_MP_14_NRGN.pdf", width=6*.9, height=5*.9, dpi = 700, units = "cm")
ggsave(plot = p3, path="02_st-xenium/figures", filename= "GEP4_gsea_plot_MP_6_MES.pdf", width=6*.9, height=5*.9, dpi = 700, units = "cm")




##  --------GEP heatmap -----------

norm_d <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/data/filtered_norm_xenium_counts.csv", row.names = 1)
d_meta <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/data/filtered_norm_xenium_counts_obs_metadata.csv", row.names = 1)
topgenes <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/xenium_cNMF_k4_topgenes.csv", row.names = 1)

topgenes <- as.list(topgenes)
names(topgenes) <- paste0("GEP", 1:4)
topgenes <- reshape2::melt(topgenes)
colnames(topgenes) <- c("gene", "GEP")

norm_d <- t(norm_d)[topgenes$gene,]
norm_d <- t(scale(t(norm_d)))

d_meta <- d_meta %>%
  mutate(GEP1_assignement = as.factor(if_else(grepl("1", GEP_assignment), 1, 0)),
         GEP2_assignement = as.factor(if_else(grepl("2", GEP_assignment), 1, 0)),
         GEP3_assignement = as.factor(if_else(grepl("3", GEP_assignment), 1, 0)),
         GEP4_assignement = as.factor(if_else(grepl("4", GEP_assignment), 1, 0)))


GEP1_colors <- c(
  "0" = "white",
  "1" = "#D73027"   # warm red
)

GEP2_colors <- c(
  "0" = "white",
  "1" = "#4575B4"   # cool blue
)

GEP3_colors <- c(
  "0" = "white",
  "1" = "#1A9850"   # green
)

GEP4_colors <- c(
  "0" = "white",
  "1" = "#762A83"   # purple
)

top_ha <- HeatmapAnnotation(
  GEP1 = d_meta$GEP1_assignement,
  GEP2 = d_meta$GEP2_assignement,
  GEP3 = d_meta$GEP3_assignement,
  GEP4 = d_meta$GEP4_assignement,
  #batch = d_meta$batch,
  col = list(GEP1 = GEP1_colors, GEP2 = GEP2_colors, GEP3 = GEP3_colors, GEP4 = GEP4_colors),
  annotation_height = unit(6, "mm")
)

# Genes you want to annotate
genes_to_mark <- c(
  # Stress / HSP / ISR program
  "HSPA6", "HSPA9", "HSPB8", "DNAJA1", "ATF4",
  
  # Mesenchymal / ECM–remodeling program
  "LOX", "MMP14", "TNC", "TGFB3", "CD44", "PDPN",
  
  # Inflammatory / cytokine–EGFR program
  "IL1B", "IL6", "CSF3", "FOSL1", "EGFR", "VEGFA",
  
  # Regulatory / splicing–epigenetic program
  "SF3B1", "SRSF1", "EZH2", "DOT1L", "HMGA2"
)

# Convert gene names to row indices in norm_d
mark_idx <- match(genes_to_mark, rownames(norm_d))

# Remove genes that were not found
valid <- !is.na(mark_idx)
mark_idx <- mark_idx[valid]
mark_labels <- genes_to_mark[valid]

# Create row annotation with anno_mark
row_ha = rowAnnotation(
  gene_mark = anno_mark(
    at = mark_idx,
    labels = mark_labels,
    labels_gp = gpar(fontsize = 8)
  )
)

pdf("02_st-xenium/figures/gep_heatmap.pdf", width = 5, height = 6)
Heatmap(
  norm_d,
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  name = "z",
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  cluster_row_slices = FALSE,
  top_annotation = top_ha,
  row_split = factor(topgenes$GEP),
  show_column_dend = FALSE,
  border = TRUE,
  right_annotation = row_ha   # <-- add mark labels on the right
)
dev.off()

## -------- GEP assignments across samples ---------

library(ggpubr)

df <- d_meta
df$batch <- gsub("table-", "", df$batch)
df$batch <- factor(df$batch)
levels(df$batch) <- c("KSC|GBM", "KSC|+M2", "KSC|+MΦ", "KSC|+Mono", "uG|GBM", "uG|+M2", "uG|+MΦ", "uG|+Mono")

p1 <- df %>% dplyr::select(starts_with("Usage_"), batch) %>% rownames_to_column() %>%
  pivot_longer(starts_with("Usage_"), names_to = "GEP", values_to = "Score") %>%
  mutate(GEP = factor(GEP, levels = paste0("Usage_", 1:4), labels = paste0("GEP", 1:4))) %>%
  ggboxplot(x = "batch", y = "Score") +
  facet_grid(. ~ GEP) +
  labs(x = NULL, y = "GEP usage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        strip.background = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

cell_counts <- df %>% group_by(batch) %>% summarise(Cell_Count = n())

p2 <- df %>% dplyr::select(ends_with("_assignement"), batch) %>% rownames_to_column() %>%
  pivot_longer(ends_with("_assignement"), names_to = "GEP", values_to = "Assignment") %>%
  group_by(batch, GEP, Assignment) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  left_join(cell_counts, "batch") %>%
  mutate(Prop = Count / Cell_Count) %>%
  filter(Assignment == 1) %>%
  mutate(GEP = gsub("_assignement", "", GEP)) %>%
  ggbarplot(x = "batch", y = "Prop", fill = "black") +
  facet_grid(. ~ GEP) +
  labs(x = NULL, y = "Frac. of cells expressing GEP") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        strip.background = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

p <- p1 / p2

ggsave(
  plot = p,
  filename = "GEP_exp_across_batches.svg",
  path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures/",
  device = svglite::svglite,
  width=15, height=15, units = "cm"
)

# Number of active GEPs

no_active_geps <- d_meta %>%
  dplyr::select(ends_with("_assignement"), batch) %>%
  dplyr::mutate(
    dplyr::across(
      ends_with("_assignement"),
      ~ as.numeric(.x) - 1
    ),
    GEP_assignments_sum = rowSums(
      dplyr::across(ends_with("_assignement")),
      na.rm = TRUE
    )
  )

df_plot <- no_active_geps %>%
  group_by(batch, GEP_assignments_sum) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(Freq = Count / sum(Count)) %>%
  ungroup() %>%
  mutate(GEP_assignments_sum = factor(GEP_assignments_sum), batch = gsub("table-", "", batch))

p <- ggplot(df_plot, aes(x = batch, y = Freq, fill = GEP_assignments_sum)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Fraction of cells", fill = "# active GEPs") +
  scale_fill_brewer(palette = "Reds") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

ggsave(plot = p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures/", filename= "active_GEP_numbers_across_batches.pdf", width=12, height=6, dpi = 700, units = "cm")

no_active_geps_summary <- no_active_geps %>%
  group_by(GEP_assignments_sum) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(Pct = round(Count / sum(Count) * 100)) %>%
  ungroup() %>%
  mutate(GEP_assignments_sum = factor(GEP_assignments_sum))

# A tibble: 4 × 3
#GEP_assignments_sum Count    Freq
#<fct>               <int>   <dbl>
#  1 0                    1009 0.193  
#2 1                    3101 0.592  
#3 2                    1107 0.211  
#4 3                      22 0.00420

##  -------- Quadrant plots -----------

### ---------Assignements ------------

usage_d <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/xenium_cNMF_k4_usage_norm.csv", row.names = 1)

M <- quadrant_matrix(usage_d)

colnames(M) <- c("Quadrant_X", "Quadrant_Y")

d_meta <- cbind(d_meta, M)

gep_col_vec <- setNames(c("#D73027", "#4575B4", "#1A9850", "#762A83"), paste0("GEP", 1:4))

p_list <- lapply(1:4, function(i) {
  x <- paste0("GEP", i, "_assignement")
  col <- gep_col_vec[i]
  ggplot(d_meta, aes(x = Quadrant_X, 
                     y = Quadrant_Y, 
                     color = as.factor(!!sym(x)))) +
    geom_point(alpha = 0.7, size = 0.2) +
    facet_wrap(~ batch, ncol = 8) +
    labs(
      x = NULL,
      y = NULL,
      color = paste0("GEP", i)
    ) +
    theme_bw() +
    scale_color_manual(values = c(
      "0" = "grey",
      "1" = col
    )) +
    theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
})
  
p <- wrap_plots(p_list, ncol = 1)
ggsave(p, path="02_st-xenium/figures", filename= "gep_quadrant_plot_w_assignments.pdf", width=10, height=6, units = "in", dpi = 700)

### ------Scatterpie-----------

library(scatterpie)
gep_col_vec <- setNames(c("#D73027", "#4575B4", "#1A9850", "#762A83"), paste0("Usage_", 1:4))




d_meta$Perturbation <- str_split(str_split(d_meta$batch, "-", Inf, T)[,2], "_", Inf, T)[,1]
d_meta$Composition <- str_split(str_split(d_meta$batch, "-", Inf, T)[,2], "_", Inf, T)[,2]

p <- ggplot() +
  geom_scatterpie(
    aes(x = Quadrant_X, y = Quadrant_Y),
    data = d_meta,
    cols = paste0("Usage_", 1:4),
    pie_scale = 0.5,
    colour = NA
  ) +
  facet_grid(Perturbation ~ Composition) +
  scale_fill_manual(values = gep_col_vec) +
  labs(
    x = NULL,
    y = NULL,
    fill = "GEP"
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(p, path="02_st-xenium/figures", filename= "gep_quadrant_plot_w_pies.pdf", width=5, height=2.25, units = "in", dpi = 700)


p <- ggplot() +
  geom_scatterpie(
    aes(x = Quadrant_X, y = Quadrant_Y),
    data = d_meta,
    cols = paste0("Usage_", 1:4),
    pie_scale = 0.3,
    colour = NA
  ) +
  scale_fill_manual(values = gep_col_vec) +
  labs(
    x = NULL,
    y = NULL,
    fill = "GEP"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(p, path="02_st-xenium/figures", filename= "gep_quadrant_plot_w_pie.pdf", width=3.25, height=2.25, units = "in", dpi = 700)

## ------Scores correlation---------

library(ComplexHeatmap)
library(circlize)

scores_df <- d_meta %>% dplyr::select(ends_with("_score"))

cor_mat <- cor(scores_df, method = "spearman")
colnames(cor_mat) <- paste0("GEP", 1:4)
rownames(cor_mat) <- paste0("GEP", 1:4)

col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("#2166AC", "white", "#B2182B")
)

ht <- Heatmap(
  cor_mat,
  col = col_fun,
  name = "rho",
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp    = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- cor_mat[i, j]
    if (is.na(v)) return()
    
    # choose text color based on how dark the tile is (simple heuristic)
    txt_col <- if (abs(v) > 0.5) "white" else "black"
    
    grid.text(
      sprintf("%.2f", v),
      x, y,
      gp = gpar(fontsize = 8, col = txt_col)
    )
  }
)

pdf("02_st-xenium/figures/GEP_cor_heatmap.pdf", width = 3*0.9, height = 2.4*0.9)
draw(ht)
dev.off()


## -------GEP frequency analysis-------

gep_assignments <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/xenium_cNMF_k4_gep_assignments.csv", row.names = 1)

gep_assignments <- gep_assignments %>%
  mutate(Composition_lvl1 = if_else(grepl("U87", batch), "U87", "U87_THP-1"),
         Composition_lvl2 = str_split(batch, "_", Inf, T)[,2],
           Perturbation = if_else(grepl("uG", batch), "uG", "KSC"))

geps <- as.character(1:4)

# crossing creates a cartesian product of two dataframes, i.e. it pairs each row of geps with each row of gep_assignments
gep_freqs <- tidyr::crossing(GEP = geps, gep_assignments) %>%
  mutate(
    is_gep = str_detect(GEP_assignment, GEP)
  ) %>%
  group_by(Composition_lvl1, Composition_lvl2, Perturbation, batch, GEP, is_gep) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(freq = count / sum(count)) %>%
  drop_na() %>%                 # removes GEPs with no matches in a batch
  filter(is_gep) %>%            # keep only TRUE
  dplyr::select(-is_gep) %>%
  mutate(GEP = paste0("GEP ", GEP))

gep_col_vec <- setNames(c("#D73027", "#4575B4", "#1A9850", "#762A83"), paste0("GEP ", 1:4))
gep_freqs$freq <- gep_freqs$freq*100

p1 <- ggplot(gep_freqs, aes(
  x = Perturbation,
  y = freq,
  group = GEP,
  color = GEP
)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  facet_grid(~ Composition_lvl2) +
  theme_cowplot(14) +
  theme(
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0.8, "lines")
  ) +
  labs(
    x = NULL,
    y = "% expressing",
    color = "Group"
  ) + 
  scale_color_manual(values = gep_col_vec)

ggsave(p1, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures", filename= "GEP_frequency.pdf", width=6*1.2, height=2*1.2, units = "in", dpi = 700)


gep_freqs %>%
  ungroup() %>%
  dplyr::select(Composition_lvl2, Perturbation, freq, GEP) %>%
  pivot_wider(names_from = "Perturbation", values_from = "freq") %>%
  mutate(FC = uG / KSC, LFC = log2(FC)) %>%
  filter(abs(LFC) > 1)

# U87: GEP4 goes up. According to GO GSEA of GEP 4, "negative regulation of gene expression, epigenetic" is enriched, which is consistent with mRNA-seq. However, GEP 2 (MES) does not go down.