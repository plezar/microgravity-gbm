library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

ref_object <- readRDS("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/data/U87_THP_combined_ref.rds")
ref_object$cell_type <- str_split(ref_object$orig.ident, "_", Inf, T)[,1]

p <- DimPlot(ref_object, group.by = "cell_type", pt.size = 0.1) +
  ggtitle(NULL) +
  scale_color_brewer(palette = "Set2") +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(plot = p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures/", filename= "scref_umap.pdf", width=9*1.2, height=7*1.2, dpi = 700, units = "cm")

Idents(ref_object) <- "cell_type"

ref_object <- JoinLayers(ref_object)
all_markers <- FindAllMarkers(ref_object, only.pos = T)


# ---- Plot canonical monocyte markers ---

## ---- scRNAseq ref ------

monocyte_markers <- c(
  "CD14","FCN1","CCR2","CSF1R","FCGR1A","FCGR3A","ITGAM",
  "ITGB2","SPI1","IRF8","SIRPA","CLEC12A","BST1"
)

# choose grouping (matches what DotPlot would use by default: Idents(ref_object))
group_var <- "ident"  # internal name we'll use in the code

# --- extract expression for just these genes (works for RNA assay by default) ---
# If you want a specific assay/slot, set these:
assay_use <- DefaultAssay(ref_object)  # e.g., "RNA"
slot_use  <- "data"                   # "data"=log-normalized, "counts"=raw counts

expr <- FetchData(
  object = ref_object,
  vars   = monocyte_markers,
  slot   = slot_use,
  assay  = assay_use
)

meta <- data.frame(ident = Idents(ref_object))
df <- cbind(meta, expr)

# pct expressed per group (expression > 0 in the chosen slot)
pct_df <- df %>%
  pivot_longer(cols = all_of(monocyte_markers),
               names_to = "gene", values_to = "expr") %>%
  group_by(ident, gene) %>%
  summarise(pct_exp = 100 * mean(expr > 0, na.rm = TRUE), .groups = "drop")

# keep the same gene order as your vector
pct_df$gene <- factor(pct_df$gene, levels = rev(monocyte_markers))

# plot: dot size = pct expressed
p1 <- ggplot(pct_df, aes(x = ident, y = gene, size = pct_exp, color = ident)) +
  geom_point(
    shape = 21,          # hollow circle
    stroke = 0.8,        # outline thickness
    fill = NA           # no fill
  ) +
  scale_color_manual(values = c("THP1" = "#66C2A5",
                     "U87"  = "#FC8D62")) +
  scale_size_continuous(limits = c(0, 100), range = c(0.1, 6)) +
  labs(x = NULL, y = NULL, size = "% expressed") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.line = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


## ---- Xenium ------

expr <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/data/filtered_norm_xenium_counts.csv", row.names = 1)
meta <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/data/filtered_norm_xenium_counts_obs_metadata.csv", row.names = 1)
df <- cbind(meta, expr)

# pct expressed per group (expression > 0 in the chosen slot)
pct_df <- df %>%
  pivot_longer(cols = all_of(monocyte_markers),
               names_to = "gene", values_to = "expr") %>%
  group_by(batch, gene) %>%
  summarise(pct_exp = 100 * mean(expr > 0, na.rm = TRUE), .groups = "drop")

# keep the same gene order as your vector
pct_df$gene <- factor(pct_df$gene, levels = rev(monocyte_markers))

pct_df$batch <- gsub("table-", "", pct_df$batch)

# plot: dot size = pct expressed
p2 <- ggplot(pct_df, aes(x = batch, y = gene, size = pct_exp)) +
  geom_point(
    shape = 21,          # hollow circle
    stroke = 0.8,        # outline thickness
    fill = NA           # no fill
  ) +
  scale_size_continuous(limits = c(0, 100), range = c(0.1, 6)) +
  labs(x = NULL, y = NULL, size = "% expressed") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

p <- p1 + p2 + plot_layout(widths = c(1, 4))

ggsave(plot = p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures/", filename= "scref_monocyte_marker_expression.pdf", width=12, height=10, dpi = 700, units = "cm")







myeloid_markers_core <- c(
  # Lineage / identity
  "TYROBP",   # DAP12
  "SPI1",     # PU.1
  "FCER1G",
  "AIF1",     # IBA1
  "GMFG",
  
  # Surface / leukocyte markers
  "PTPRC",    # CD45
  "ITGB2",    # CD18
  "MS4A7",
  "ANPEP",    # CD13
  "CD53",
  
  # Innate immune / phagocytic machinery
  "CTSZ",
  "VAMP8",
  "NCF4",
  "PYCARD",
  
  # Macrophage metabolism / activation
  "APOE",
  "PLIN2"
)

gbm_mes_markers <- c(
  # Core MES GBM identity
  "ALDH1A3",   # hallmark MES GBM marker
  "SERPINE1",  # MES / hypoxia / invasion
  "TNC",       # ECM, MES GBM
  "SPARC",     # invasion / ECM remodeling
  "LOXL2",     # ECM crosslinking, MES state
  
  # MES-associated signaling / structure
  "CALD1",     # cytoskeleton, mesenchymal
  "AKAP12",    # mesenchymal-like GBM
  "TFPI2",     # invasion / ECM regulation
  "SRPX",      # mesenchymal / astrocytic GBM
  "TM4SF1"     # tumor cell motility
)

all_markers <- unique(c(gbm_mes_markers, myeloid_markers_core))

E <- GetAssayData(ref_object, slot = "data")
E <- as.matrix(E[all_markers,])
E <- t(scale(t(E)))



# make sure ordering matches the heatmap columns
cell_type <- ref_object$cell_type[colnames(E)]

pal <- RColorBrewer::brewer.pal(name = "Set2", n = 3)
# define colors (adjust as needed)
cell_type_cols <- c(
  "THP1" = "#66C2A5",
  "U87"  = "#FC8D62"
)



top_ha <- HeatmapAnnotation(
  CellType = cell_type,
  col = list(CellType = cell_type_cols),
  show_annotation_name = TRUE
)

pdf("02_st-xenium/figures/RCTD_ref_markers.pdf", width = 8*.6, height = 5*.6)
Heatmap(
  E,
  col = colorRamp2(
    c(-3, 0, 3),
    c("#2166AC", "white", "#B2182B")
  ),
  top_annotation = top_ha,
  show_column_dend = FALSE,
  border = T,
  row_names_gp = grid::gpar(fontsize = 6),
  show_row_dend = FALSE,
  show_column_names = FALSE,
  name = "z"
)
dev.off()

