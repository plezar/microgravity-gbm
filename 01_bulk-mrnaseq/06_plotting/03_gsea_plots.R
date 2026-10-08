library(dplyr)
library(tibble)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(enrichplot)
library(ggplot2)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
source("utils/plot_gsea_annot.R")
source("utils/run_gsea_cp.R")


leading_edge_list <- function(gsea_cp, gsets = NULL, split_pat = "/") {
  res <- as.data.frame(gsea_cp@result)  # or as.data.frame(gsea_cp)
  
  if (!is.null(gsets)) {
    res <- res %>% filter(ID %in% gsets)
  }
  
  # core_enrichment is a single string like "GENE1/GENE2/GENE3"
  le <- res %>%
    transmute(
      ID,
      leading_edge = strsplit(core_enrichment, split = split_pat, fixed = TRUE)
    )
  
  # named list: names are pathway IDs, values are character vectors of genes
  out <- le$leading_edge
  names(out) <- le$ID
  out
}

make_le_heatmap <- function(gsea_cp,
                            dds,
                            de_stats,
                            gsets = c("MP_1_RP", "MP_5_Hypoxia", "MP_15_Stress2"),
                            cell_group = "U87",
                            gravity_levels = c("KSC", "uG"),
                            anno_size_mm = 2.5,
                            gravity_colors = c(KSC = "grey", uG = "black"),
                            lfc_limits = c(-3, 3),
                            lfc_palette = "PRGn",
                            row_name_fontsize = 7,
                            pdf_file = NULL,
                            pdf_width = 3.5,
                            pdf_height = 4) {
  stopifnot(!missing(gsea_cp), !missing(dds), !missing(de_stats))
  
  # deps
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("stringr", quietly = TRUE)
  requireNamespace("ComplexHeatmap", quietly = TRUE)
  requireNamespace("circlize", quietly = TRUE)
  requireNamespace("RColorBrewer", quietly = TRUE)
  requireNamespace("grid", quietly = TRUE)
  requireNamespace("DESeq2", quietly = TRUE)
  
  # --- leading edge genes + row split factor (MP_1, MP_5, MP_15) ---
  le_genes <- leading_edge_list(gsea_cp, gsets)
  le_all <- unlist(le_genes, use.names = FALSE)
  
  group_factor <- rep(names(le_genes), lengths(le_genes))
  group_factor <- stringr::str_split_fixed(group_factor, "_", 3)[, 1:2]
  group_factor <- apply(group_factor, 1, paste0, collapse = "_")
  group_factor <- factor(group_factor, levels = unique(group_factor))
  
  # --- expression matrix (VST) ---
  vsd <- DESeq2::vst(dds)
  vsd_mat <- SummarizedExperiment::assay(vsd)
  
  # ensure genes exist
  le_all <- intersect(le_all, rownames(vsd_mat))
  if (length(le_all) == 0) stop("None of the leading-edge genes are present in vst(dds).")
  
  # subset columns by cell group
  cd <- SummarizedExperiment::colData(dds)
  if (!("Group_cell" %in% colnames(cd))) stop("colData(dds)$Group_cell not found.")
  if (!("Group_gravity" %in% colnames(cd))) stop("colData(dds)$Group_gravity not found.")
  
  keep_cols <- cd$Group_cell == cell_group
  if (!any(keep_cols)) stop("No samples found with Group_cell == ", cell_group)
  
  # z-score by gene (row)
  M <- t(scale(t(vsd_mat[le_all, , drop = FALSE])))
  M_sub <- M[, keep_cols, drop = FALSE]
  
  # --- top annotation ---
  grp <- factor(cd$Group_gravity[keep_cols], levels = gravity_levels)
  
  ha_col <- ComplexHeatmap::HeatmapAnnotation(
    Group_gravity = grp,
    annotation_name_side = "left",
    simple_anno_size = grid::unit(anno_size_mm, "mm"),
    annotation_name_gp = grid::gpar(fontsize = 0),
    col = list(Group_gravity = gravity_colors)
  )
  
  # --- LFC side heatmap ---
  if (!("log2FoldChange" %in% colnames(de_stats))) {
    stop("de_stats must contain a 'log2FoldChange' column.")
  }
  missing_lfc <- setdiff(le_all, rownames(de_stats))
  if (length(missing_lfc) > 0) {
    warning("Some LE genes missing in de_stats and will be dropped from LFC heatmap: ",
            paste(missing_lfc, collapse = ", "))
  }
  lfc_genes <- intersect(le_all, rownames(de_stats))
  
  lfc_mat <- as.matrix(de_stats[lfc_genes, "log2FoldChange", drop = FALSE])
  col_fun <- circlize::colorRamp2(
    seq(lfc_limits[1], lfc_limits[2], length.out = 9),
    RColorBrewer::brewer.pal(9, lfc_palette)
  )
  
  lfc_ht <- ComplexHeatmap::Heatmap(
    lfc_mat,
    rect_gp = grid::gpar(col = "white", lwd = 2),
    col = col_fun,
    row_names_gp = grid::gpar(fontsize = row_name_fontsize),
    show_column_names = FALSE,
    name = "LFC"
  )
  
  # Align row order between main and LFC heatmaps by using the same row order input
  # (main uses le_all; lfc uses lfc_genes, so subset main to lfc_genes to keep rows identical)
  # If you prefer to keep all rows in main even if LFC is missing, remove this block.
  common_rows <- intersect(rownames(M_sub), rownames(lfc_mat))
  M_sub <- M_sub[common_rows, , drop = FALSE]
  lfc_mat <- lfc_mat[common_rows, , drop = FALSE]
  # rebuild lfc_ht after subsetting
  lfc_ht <- ComplexHeatmap::Heatmap(
    lfc_mat,
    rect_gp = grid::gpar(col = "white", lwd = 2),
    col = col_fun,
    row_names_gp = grid::gpar(fontsize = row_name_fontsize),
    show_column_names = FALSE,
    name = "LFC"
  )
  
  # subset group_factor to common rows
  # (group_factor corresponds to le_all order; rebuild for common_rows)
  idx_map <- match(common_rows, le_all)
  group_factor2 <- group_factor[idx_map]
  group_factor2 <- factor(group_factor2, levels = unique(group_factor2))
  
  main_ht <- ComplexHeatmap::Heatmap(
    M_sub,
    top_annotation = ha_col,
    row_split = group_factor2,
    cluster_columns = FALSE,
    name = "z",
    border = TRUE,
    show_row_dend = FALSE,
    show_column_names = FALSE,
    row_names_gp = grid::gpar(fontsize = row_name_fontsize)
  )
  
  ht <- main_ht + lfc_ht
  
  if (!is.null(pdf_file)) {
    grDevices::pdf(pdf_file, width = pdf_width, height = pdf_height)
    ComplexHeatmap::draw(ht)
    grDevices::dev.off()
  }
  
  invisible(list(
    ht = ht,
    main_ht = main_ht,
    lfc_ht = lfc_ht,
    le_genes = le_genes,
    group_factor = group_factor2,
    M_sub = M_sub,
    lfc_mat = lfc_mat
  ))
}

# Example:
# res <- make_le_heatmap(out$gsea_cp, dds, de_stats,
#                        gsets = c("MP_1_RP","MP_5_Hypoxia","MP_15_Stress2"),
#                        cell_group = "U87",
#                        pdf_file = "01_bulk-mrnaseq/results/figures/U87_LE_heatmap.pdf")
# ComplexHeatmap::draw(res$ht)

#################
# === U87 ===
#################

# === Read MP sets ===
sig <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/data/Nomura_2025_signatures.csv")
sig <- sig[,grepl("^MP", colnames(sig))]

# === DE data ====
de_d <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87.csv", row.names = 1)

# === GSEA ===
out <- run_gsea_cp(sig, de_d)

# === Plot ===
as.data.frame(out$gsea_cp)

p1 <- plot_gsea_annot(out$gsea_cp, "MP_1_RP", color = "#7b0f0f", label_size=4, label_alpha=0.5)
p2 <- plot_gsea_annot(out$gsea_cp, "MP_5_Hypoxia", color = "#0b3c6f", label_size=4, label_alpha=0.5)
p3 <- plot_gsea_annot(out$gsea_cp, "MP_6_MES", color = "#0b3c6f", label_size=4, label_alpha=0.5)

ggsave(plot = p1, path="01_bulk-mrnaseq/results/figures", filename= "U87_MP_1_RP_gsea_plot.pdf", width=7, height=5, dpi = 700, units = "cm")
ggsave(plot = p2, path="01_bulk-mrnaseq/results/figures", filename= "U87_MP_5_Hypoxia_gsea_plot.pdf", width=7, height=5, dpi = 700, units = "cm")
ggsave(plot = p3, path="01_bulk-mrnaseq/results/figures", filename= "U87_MP_15_Stress2_gsea_plot.pdf", width=7, height=5, dpi = 700, units = "cm")

# === Leading edge analysis ===

dds <- readRDS("01_bulk-mrnaseq/results/dds.rds")
de_stats <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87.csv", row.names = 1)
res <- make_le_heatmap(out$gsea_cp, dds, de_stats,
                       gsets = c("MP_1_RP","MP_5_Hypoxia","MP_6_MES"),
                       cell_group = "U87",
                       row_name_fontsize = 5,
                       pdf_width = 3.5,
                       pdf_file = "01_bulk-mrnaseq/results/figures/U87_LE_heatmap.pdf")

#pdf("01_bulk-mrnaseq/results/figures/U87_LE_heatmap.pdf", width = 3.5, height = 4)
#draw(main_ht + lfc_ht)
#dev.off()

#################
# === Neftel 2019 ===
#################

# === Read MP sets ===
sig_neftel <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/data/Neftel_2019_signatures.csv")

# === DE data ====
de_d <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87.csv", row.names = 1)

# === GSEA ===
out <- run_gsea_cp(sig_neftel, de_d)

# === Plot ===
as.data.frame(out$gsea_cp)

p1 <- plot_gsea_annot(out$gsea_cp, "MES2", color = "#0b3c6f", label_size=4, label_alpha=0.5)
ggsave(plot = p1, path="01_bulk-mrnaseq/results/figures", filename= "U87_MES2_gsea_plot.pdf", width=7, height=5, dpi = 700, units = "cm")


###################
# === U87 THP1 ===
###################

# === DE data ====
de_d <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87THP1.csv", row.names = 1)

# === GSEA ===
out <- run_gsea_cp(sig, de_d)

# === Plot ===
as.data.frame(out$gsea_cp)

p1 <- plot_gsea_annot(out$gsea_cp, "MP_1_RP", color = "grey75", label_size=4, label_alpha=0.5)
p2 <- plot_gsea_annot(out$gsea_cp, "MP_5_Hypoxia", color = "grey75", label_size=4, label_alpha=0.5)
p3 <- plot_gsea_annot(out$gsea_cp, "MP_15_Stress2", color = "#0b3c6f", label_size=4, label_alpha=0.5)

ggsave(plot = p1, path="01_bulk-mrnaseq/results/figures", filename= "U87THP1_MP_1_RP_gsea_plot.pdf", width=7, height=5, dpi = 700, units = "cm")
ggsave(plot = p2, path="01_bulk-mrnaseq/results/figures", filename= "U87THP1_MP_5_Hypoxia_gsea_plot.pdf", width=7, height=5, dpi = 700, units = "cm")
ggsave(plot = p3, path="01_bulk-mrnaseq/results/figures", filename= "U87THP1_MP_15_Stress2_gsea_plot.pdf", width=7, height=5, dpi = 700, units = "cm")


#################
# === Neftel 2019 ===
#################

# === GSEA ===
out <- run_gsea_cp(sig_neftel, de_d)

# === Plot ===
as.data.frame(out$gsea_cp)

p1 <- plot_gsea_annot(out$gsea_cp, "MES2", color = "#0b3c6f", label_size=4, label_alpha=0.5)
ggsave(plot = p1, path="01_bulk-mrnaseq/results/figures", filename= "U87THP1_MES2_2019_gsea_plot.pdf", width=7, height=5, dpi = 700, units = "cm")

# === HALLMARK ===

gsea_res <- readRDS("01_bulk-mrnaseq/results/enrichment/deseq2_uG_vs_KSC_U87THP1/GSEA/Hallmark_GSEA.rds")

## === Plot ===
#as.data.frame(gsea_res@result) %>% View()
gsets <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_COMPLEMENT")

p1 <- plot_gsea_annot(gsea_res, gsets[1], color = "#7b0f0f", label_size=4, label_alpha=0.5, wrap_width = 18)
p2 <- plot_gsea_annot(gsea_res, gsets[2], color = "#7b0f0f", label_size=4, label_alpha=0.5, wrap_width = 18)
p3 <- plot_gsea_annot(gsea_res, gsets[3], color = "#7b0f0f", label_size=4, label_alpha=0.5, wrap_width = 18)

ggsave(plot = p1, path="01_bulk-mrnaseq/results/figures", filename= "U87THP1_TNFA_gsea_plot.pdf", width=7, height=5, dpi = 700, units = "cm")
ggsave(plot = p2, path="01_bulk-mrnaseq/results/figures", filename= "U87THP1_KRAS_gsea_plot.pdf", width=7, height=5, dpi = 700, units = "cm")
ggsave(plot = p3, path="01_bulk-mrnaseq/results/figures", filename= "U87THP1_Compl_gsea_plot.pdf", width=7, height=5, dpi = 700, units = "cm")


## === Graph Plot ===
de_stats <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87THP1.csv", row.names = 1)
gsea_res <- clusterProfiler::setReadable(gsea_res, "org.Hs.eg.db", keyType = "ENTREZID")

fcs <- setNames(de_stats$log2FoldChange, rownames(de_stats))
genes_to_label <- c(
  "CCL20","IL1B","GFPT2","BIRC3",
  "ID2","PLAUR","IL6","SERPINE1"
)

p1 <- cnetplot(gsea_res, foldChange=fcs, node_label="all", color_category = "black", color_edge = "grey75") +
  scale_color_gradient2(name='LFC', low='darkblue', high='firebrick')
p1$layers[[4]] <- NULL
p2 <- p1 + ggraph::geom_node_text(data=p1$data[which(p1$data$name %in% c(p1$data$name[1:4], unname(genes_to_label)) ),], #c(p$data$name[1:3],unname(unlist(y_short))) ),],
                         aes(x=x, y=y, label=name), color='black', bg.color="white",
                         segment.size=0.5, repel=TRUE, size=3, min.segment.length = 0, max.overlaps = 20)

ggsave(plot = p2, path="01_bulk-mrnaseq/results/figures", filename= "U87THP1_gsea_cnetplot.pdf", width=12, height=7, dpi = 700, units = "cm")


### === Testing for overlap ====

hallmarks <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
               "HALLMARK_INFLAMMATORY_RESPONSE",
               "HALLMARK_KRAS_SIGNALING_UP",
               "HALLMARK_COMPLEMENT")

# --- 1) Leading-edge sets from your GSEA result ---
le_df <- as.data.frame(gsea_res) %>%
  filter(ID %in% hallmarks) %>%
  dplyr::select(ID, core_enrichment) %>%
  mutate(le_genes = str_split(core_enrichment, "/"))

LE <- setNames(le_df$le_genes, le_df$ID)

# --- 2) Full Hallmark gene sets from MSigDB (match species to your data) ---
# Choose ONE:
species <- "Homo sapiens"  # or "Mus musculus"

library(msigdbr)
msig_h <- msigdbr(species = species, category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

FULL <- msig_h %>%
  filter(gs_name %in% hallmarks) %>%
  group_by(gs_name) %>%
  summarise(full_genes = list(unique(gene_symbol)), .groups = "drop") %>%
  { setNames(.$full_genes, .$gs_name) }

# --- sanity check: all sets present ---
missing_full <- setdiff(hallmarks, names(FULL))
if (length(missing_full) > 0) {
  stop("Missing Hallmark sets in msigdbr: ", paste(missing_full, collapse = ", "))
}

# --- 3) Pairwise test: is LE overlap rate enriched beyond baseline Hallmark overlap? ---
pairs <- combn(hallmarks, 2, simplify = FALSE)

overlap_enrichment_test <- function(A, B) {
  leA <- unique(LE[[A]]); leB <- unique(LE[[B]])
  sA  <- unique(FULL[[A]]); sB <- unique(FULL[[B]])
  
  # Leading-edge overlap counts
  a_le <- length(intersect(leA, leB))
  u_le <- length(union(leA, leB))
  b_le <- u_le - a_le
  
  # Full-set overlap counts (baseline)
  a_full <- length(intersect(sA, sB))
  u_full <- length(union(sA, sB))
  b_full <- u_full - a_full
  
  # 2x2 comparing "shared vs not shared" in LE vs FULL
  mat <- matrix(c(a_le, b_le,
                  a_full, b_full),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("LeadingEdge", "FullSet"),
                                c("Shared", "NotShared")))
  
  ft <- fisher.test(mat)  # LE overlap rate > baseline overlap rate
  
  # Odds ratio: (a_le/b_le) / (a_full/b_full)
  # fisher.test$estimate corresponds to this table's OR
  tibble(
    term_A = A,
    term_B = B,
    shared_LE = a_le,
    union_LE = u_le,
    shared_full = a_full,
    union_full = u_full,
    overlap_frac_LE = ifelse(u_le > 0, a_le / u_le, NA_real_),
    overlap_frac_full = ifelse(u_full > 0, a_full / u_full, NA_real_),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value
  )
}

res2 <- map_dfr(pairs, \(p) overlap_enrichment_test(p[1], p[2])) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj, p_value)

res2

library(dplyr)
library(ggplot2)
library(scales)
library(stringr)

# res2 = your results tibble (term_A, term_B, odds_ratio, p_adj, overlap_frac_LE, ...)

plot_df <- res2 %>%
  mutate(pair = str_replace_all(paste0(term_A, "  ×  ", term_B), "HALLMARK_", "")) %>%
  mutate(pair = factor(pair, levels = pair[order(odds_ratio, decreasing = TRUE)])) %>%
  mutate(sig = case_when(
    p_adj < 0.01 ~ "**",
    p_adj < 0.05 ~ "*",
    p_adj < 0.10 ~ "·",
    TRUE ~ ""
  ))

ggplot(plot_df, aes(x = odds_ratio, y = pair)) +
  geom_segment(aes(x = 1, xend = odds_ratio, yend = pair), linewidth = 0.6) +
  geom_point(aes(size = overlap_frac_LE, color = p_adj), alpha = 0.95) +
  geom_text(aes(label = paste0("OR=", round(odds_ratio, 2), sig)),
            hjust = -0.05, size = 3.2) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.5) +
  scale_color_continuous(trans = "reverse", breaks = c(0.01, 0.05, 0.10, 0.20),
                         labels = label_number(accuracy = 0.01)) +
  scale_size_continuous(range = c(2.5, 6), labels = label_percent(accuracy = 0.1)) +
  coord_cartesian(xlim = c(0.8, max(plot_df$odds_ratio) * 1.25)) +
  labs(x = "Odds ratio: LE overlap enrichment vs baseline Hallmark overlap",
       y = NULL,
       color = "FDR (BH)",
       size = "LE overlap\nfraction") +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 9),
        legend.position = "right")


# === Leading edge analysis ===

de_stats <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87THP1.csv", row.names = 1)
res <- make_le_heatmap(out$gsea_cp, dds, de_stats,
                       gsets = c("MP_1_RP","MP_5_Hypoxia","MP_15_Stress2"),
                       cell_group = "U87THP",
                       pdf_file = "01_bulk-mrnaseq/results/figures/U87THP1_LE_heatmap.pdf")
