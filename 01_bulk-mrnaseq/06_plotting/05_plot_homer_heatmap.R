library(ComplexHeatmap)
library(tidyverse)
library(org.Hs.eg.db)
source("01_bulk-mrnaseq/util.R")
library(scales)
library(patchwork)
library(clusterProfiler)
library(enrichplot)
library(cowplot)

# ==== TEAD ======

extract_offset <- function(x) {
  # often looks like: "AGGAATG,+,0.83" or "-123(AGGA..., +, 0.7)" depending on HOMER version
  # grab first integer in the string
  as.integer(str_extract(x, "-?\\d+"))
}

mat <- read_tsv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/homer_out/DEG_TEAD_annotation_ghist.txt")
mat <- as.data.frame(mat)
rownames(mat) <- mat$Gene
mat$Gene <- NULL
mat <- as.matrix(mat)
colnames(mat) <- str_split(colnames(mat), "\\.\\.\\.", n = Inf, simplify = TRUE)[,1]
tfs <- c(rep("TEAD1",3), rep("TEAD2",3))

peak_anno <- read_tsv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/homer_out/DEG_TEAD_annotation.txt")



# === TSS distribution ====

colnames(peak_anno)[1] <- "PeakID"

motif_cols <- names(peak_anno)[str_detect(names(peak_anno), "^TEAD")]
pa_long <- peak_anno %>%
  mutate(has_tead = if_any(all_of(motif_cols), ~ !is.na(.) & . != "")) %>%
  filter(has_tead) %>%
  select(PeakID, `Gene Name`, `Distance to TSS`, all_of(motif_cols)) %>%
  pivot_longer(cols = all_of(motif_cols), names_to = "motif", values_to = "motif_hit") %>%
  filter(!is.na(motif_hit) & motif_hit != "") %>%
  mutate(motif_offset = vapply(motif_hit, extract_offset, integer(1)))

p <- ggplot(pa_long, aes(x = motif_offset)) +
  geom_density(fill = "#F6E27F", alpha = 0.3, linewidth = 0.5, adjust = 0.5) +  # >1 smoother, <1 bumpier
  geom_vline(xintercept = 0, linetype = 2) +
  theme_cowplot() +
  labs(x = "TSS dist. (bp)", y = "Density") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth =1),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave(plot = p, path="01_bulk-mrnaseq/results/figures", filename= "HOMER_TEAD_distribution.pdf", width=6*1.25, height=3*1.25, dpi = 700, units = "cm")


# === ORA ====

count_mtx <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/data/count_mtx.csv", row.names = 1)
all_universe <- mapIds(
  org.Hs.eg.db,
  keys = rownames(count_mtx),
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)

set <- peak_anno$`Entrez ID`[rowSums(is.na(peak_anno[,motif_cols]))>0]

universe <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(count_mtx),
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)

go_ora_down <- enrichGO(gene          = set,
                      universe      = universe,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)


go_ora_down <- simplify(go_ora_down)

gos_to_plot <- c("GO:0001837", "GO:0002053", "GO:0010463", "GO:0040036", "GO:2000648", "GO:0001570", "GO:0010906")

go_ora_down_select <- go_ora_down@result %>% filter(ID %in% gos_to_plot)

res_down <- prep_enrich(go_ora_down_select, n = 10, ontology = "BP")

p_down <- make_panel(
  res_down, direction = "up",
  fill_low = "#d6e6f2", fill_high = "#0b3c6f",
  #fill_low = black_low, fill_high = black_high,
  legend_pos = "bottom",
  show_y_right = FALSE
) + xlab(NULL) + scale_x_continuous(
  breaks = seq(0, 100, 20),
  labels = seq(0, 100, 20),
  expand = expansion(mult = c(0.02, 0.02)),
  limits = c(0, 100)
) + theme(legend.position = "bottom", axis.title = element_text(size = 12))

ggsave(p_down + theme(legend.position = "none"), path="01_bulk-mrnaseq/results/figures", filename= "TEAD_ORA_barplot.pdf", width=11, height=5, dpi = 700, units = "cm")
ggsave(p_down, path="01_bulk-mrnaseq/results/figures", filename= "TEAD_ORA_barplot_legend.pdf", width=11, height=8, dpi = 700, units = "cm")




# ==== SMAD ======

peak_anno <- read_tsv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/01_bulk-mrnaseq/results/homer_out/DEG_SMAD3_annotation.txt")

colnames(peak_anno)[1] <- "PeakID"

motif_cols <- names(peak_anno)[str_detect(names(peak_anno), "^Smad")]
pa_long <- peak_anno %>%
  mutate(has_tead = if_any(all_of(motif_cols), ~ !is.na(.) & . != "")) %>%
  filter(has_tead) %>%
  select(PeakID, `Gene Name`, `Distance to TSS`, all_of(motif_cols)) %>%
  pivot_longer(cols = all_of(motif_cols), names_to = "motif", values_to = "motif_hit") %>%
  filter(!is.na(motif_hit) & motif_hit != "") %>%
  mutate(motif_offset = vapply(motif_hit, extract_offset, integer(1)))

p <- ggplot(pa_long, aes(x = motif_offset)) +
  geom_density(fill = "#F6E27F", alpha = 0.3, linewidth = 0.5, adjust = 0.5) +  # >1 smoother, <1 bumpier
  geom_vline(xintercept = 0, linetype = 2) +
  theme_cowplot() +
  labs(x = "TSS dist. (bp)", y = "Density") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth =1),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave(plot = p, path="01_bulk-mrnaseq/results/figures", filename= "HOMER_SMAD_distribution.pdf", width=6*1.25, height=3*1.25, dpi = 700, units = "cm")


# === ORA ====

set <- peak_anno$`Entrez ID`[rowSums(is.na(peak_anno[,motif_cols]))>0]

universe <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(count_mtx),
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)


go_ora_up <- enrichGO(gene          = set,
                      universe      = universe,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)

res_up <- prep_enrich(go_ora_up@result, n = 2, ontology = "BP")

p_up <- make_panel(
  res_up, direction = "up",
  fill_low = "#f6d8cb", fill_high = "#7b0f0f",
  #fill_low = blue_low, fill_high = blue_high,
  legend_pos = "bottom",
  show_y_right = FALSE
) + xlab("Fold Enrichment") + scale_x_continuous(
  breaks = seq(0, 100, 20),
  labels = seq(0, 100, 20),
  expand = expansion(mult = c(0.02, 0.02)),
  limits = c(0, 100)
) + theme(legend.position = "bottom", axis.title = element_text(size = 12))

ggsave(p_up + theme(legend.position = "none"), path="01_bulk-mrnaseq/results/figures", filename= "Smad_ORA_barplot.pdf", width=11, height=3, dpi = 700, units = "cm")
ggsave(p_up, path="01_bulk-mrnaseq/results/figures", filename= "Smad_ORA_barplot_legend.pdf", width=11, height=6, dpi = 700, units = "cm")
