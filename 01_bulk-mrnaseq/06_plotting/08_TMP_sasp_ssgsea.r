suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

senescent_msc_secretome <- c(
  "ACTB", "ALB", "ACTA2", "CTGF", "YWHAE", "YWHAQ", "ANXA2P2",
  "TUBB3", "ATP5B", "DSTN", "LTF", "SOD1", "HSPA6", "SEC23A",
  "RPSA", "PSMA4", "ITIH3", "PSMA2", "MYL6", "PSMA6", "CNN3",
  "PPP1CB", "LYZ", "PPP1CA", "PSMA5", "B2M", "HBA1", "ARPC3",
  "UBE2L3", "PABPC1", "HBB", "RHOC", "MTPN", "PCBP1", "RAB11A",
  "NQO2", "AMY1A", "NUTF2", "UBE2N", "MESDC2", "HMGN2", "IGHG1",
  "RAN", "RPS13", "RPS3A", "S100A7", "SELM", "RPS25",
  "CDC42", "CYR61", "TCEB2", "PSMC2", "GREM1", "RPS12", "ATP5A1",
  "PRPS2", "SNX3", "CALML5", "EIF2S3", "RPS16", "C1orf123", "PTRF",
  "AP2S1", "VTI1B", "IGLL5", "CREG1", "APOB", "VDAC1", "PROSC",
  "RPS18", "RPS15A", "HSPB7", "CALML3", "UBE2D2", "C21orf33",
  "IGLC2", "GMFB", "NHLRC3", "STMN1", "MAP2K1", "CPNE1", "C5",
  "PAK2", "RPS8", "RAB3B", "PIN1", "TMSB10", "UBE2D1", "NENF",
  "DUSP3", "MAPRE3", "NDRG1", "ITM2B", "GNB4", "SNX6"
)

results_dir <- "01_bulk-mrnaseq/results"
figures_dir <- file.path(results_dir, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

score_name <- "TMP_SASP"
gene_set <- unique(senescent_msc_secretome)

dds <- readRDS(file.path(results_dir, "dds.rds"))
vsd <- DESeq2::vst(dds, blind = FALSE)
expr_mat <- SummarizedExperiment::assay(vsd)

genes_present <- intersect(gene_set, rownames(expr_mat))
genes_missing <- setdiff(gene_set, rownames(expr_mat))

if (length(genes_present) == 0) {
  stop("None of the TMP/SASP genes are present in vst(dds).")
}

ssgsea_param <- GSVA::ssgseaParam(expr_mat, geneSets = list(TMP_SASP = genes_present))
ssgsea_scores <- GSVA::gsva(ssgsea_param, verbose = FALSE)

sample_df <- as.data.frame(SummarizedExperiment::colData(dds)) %>%
  #rownames_to_column("Sample") %>%
  mutate(
    Group_cell = factor(Group_cell, levels = c("U87", "U87THP")),
    Group_gravity = factor(Group_gravity, levels = c("KSC", "uG"))
  )

score_df <- sample_df %>%
  mutate(
    TMP_SASP = as.numeric(ssgsea_scores[score_name, Sample]),
    score_z = as.numeric(scale(TMP_SASP))
  ) %>%
  arrange(Group_cell, Group_gravity, Sample)

comparison_df <- score_df %>%
  group_by(Group_cell) %>%
  group_modify(~{
    ksc_scores <- .x %>% filter(Group_gravity == "KSC") %>% pull(TMP_SASP)
    ug_scores <- .x %>% filter(Group_gravity == "uG") %>% pull(TMP_SASP)
    wt <- wilcox.test(ug_scores, ksc_scores, exact = FALSE)

    tibble(
      n_KSC = length(ksc_scores),
      n_uG = length(ug_scores),
      mean_KSC = mean(ksc_scores),
      mean_uG = mean(ug_scores),
      median_KSC = median(ksc_scores),
      median_uG = median(ug_scores),
      delta_mean_uG_minus_KSC = mean(ug_scores) - mean(ksc_scores),
      delta_median_uG_minus_KSC = median(ug_scores) - median(ksc_scores),
      p_value = wt$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = paste0("Wilcoxon p = ", signif(p_value, 3))
  )

gene_set_df <- tibble(
  gene_set = score_name,
  n_input = length(gene_set),
  n_present = length(genes_present),
  n_missing = length(genes_missing),
  present_genes = paste(genes_present, collapse = ";"),
  missing_genes = paste(genes_missing, collapse = ";")
)

annotation_df <- score_df %>%
  group_by(Group_cell) %>%
  summarise(
    ymin = min(TMP_SASP, na.rm = TRUE),
    ymax = max(TMP_SASP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(comparison_df, by = "Group_cell") %>%
  mutate(
    x = 1.5,
    y = ymax + pmax((ymax - ymin) * 0.14, 0.02)
  )

p <- ggplot(score_df, aes(x = Group_gravity, y = TMP_SASP, fill = Group_gravity)) +
  geom_boxplot(width = 0.28, outlier.shape = NA, color = "black", linewidth = 0.4) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 2.5,
    alpha = 0.9
  ) +
  geom_text(
    data = annotation_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  facet_wrap(~Group_cell, nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "SASP ssGSEA scores",
    x = NULL,
    y = "ssGSEA enrichment score"
  ) +
  theme_cowplot(12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.line = element_blank(),
    plot.margin = margin(5.5, 14, 5.5, 5.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

write.csv(
  gene_set_df,
  file.path(results_dir, "TMP_SASP_ssgsea_gene_set_overlap.csv"),
  row.names = FALSE
)
write.csv(
  score_df,
  file.path(results_dir, "TMP_SASP_ssgsea_scores.csv"),
  row.names = FALSE
)
write.csv(
  comparison_df,
  file.path(results_dir, "TMP_SASP_ssgsea_group_comparison.csv"),
  row.names = FALSE
)

ggsave(
  plot = p,
  path = figures_dir,
  filename = "TMP_SASP_ssgsea_scores_by_gravity.pdf",
  width = 12,
  height = 10,
  dpi = 700,
  units = "cm"
)

# Create heatmap of SASP gene expression
hm_mat <- expr_mat[genes_present, colnames(expr_mat)]
hm_mat <- hm_mat[, match(rownames(sample_df), colnames(hm_mat))]

col_annot <- HeatmapAnnotation(
  Group_cell = sample_df$Group_cell,
  Group_gravity = sample_df$Group_gravity,
  col = list(
    Group_cell = c("U87" = "#1f77b4", "U87THP" = "#ff7f0e"),
    Group_gravity = c("KSC" = "#2ca02c", "uG" = "#d62728")
  )
)

hm <- Heatmap(
  t(scale(t(hm_mat))),
  name = "VST expression",
  top_annotation = col_annot,
  show_column_names = TRUE,
  show_row_names = TRUE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  column_title = "SASP Gene Expression",
  row_title = paste0("SASP Genes (n=", nrow(hm_mat), ")"),
  column_title_side = "top"
)

pdf(
  file.path(figures_dir, "TMP_SASP_ssgsea_heatmap.pdf"),
  width = 14,
  height = 16
)
draw(hm)
dev.off()

cat(
  "TMP/SASP genes detected:", length(genes_present), "of", length(gene_set), "\n",
  "Scores saved to:", file.path(results_dir, "TMP_SASP_ssgsea_scores.csv"), "\n",
  "Comparisons saved to:", file.path(results_dir, "TMP_SASP_ssgsea_group_comparison.csv"), "\n",
  "Figure saved to:", file.path(figures_dir, "TMP_SASP_ssgsea_scores_by_gravity.pdf"), "\n"
)
