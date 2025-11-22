library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
source("02_st-xenium/util.R")

## -------Enrichment analysis-------

top_genes <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/xenium_cNMF_k4_topgenes.csv", row.names = 1)
gep_scores <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/xenium_cNMF_k4_gep_scores.csv", row.names = 1)

# res <- run_go_ora(top_genes$X1, rownames(gep_scores), species="human")
#
# res_gep1 <- run_go_gsea(setNames(gep_scores$X1, rownames(gep_scores)), species = "human", p_cutoff=1)
# res_gep2 <- run_go_gsea(setNames(gep_scores$X2, rownames(gep_scores)), species = "human", p_cutoff=1)
# res_gep3 <- run_go_gsea(setNames(gep_scores$X3, rownames(gep_scores)), species = "human", p_cutoff=1)
res_gep4 <- run_go_gsea(setNames(gep_scores$X4, rownames(gep_scores)), species = "human", p_cutoff=1)

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

pal <- brewer.pal(11, "RdBu")
pal <- rev(pal)
col_fun <- colorRamp2(
  breaks = seq(min(NES_M), max(NES_M), length.out = length(pal)),
  colors = pal
)

ht <- Heatmap(
  NES_M,
  name = "NES",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_row_dend = FALSE,
  border = T,
  # specify size here
  width = unit(2, "cm"),
  height = unit(15, "cm"),
  show_column_names = TRUE,
  cell_fun = function(j, i, x, y, w, h, fill){
    if(P_M[i, j] < 0.1) {
      grid.text('*', x, y)
    }
  }
)

pdf("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures/GEP_enrichment.pdf", width = 8, height = 8)
draw(ht, heatmap_legend_side = "top")
dev.off()

# GEP 1: Inflam.
# GEP 2: ECM.
# GEP 3: Prolif/DDR
# GEP 4: 

# res_gep1 <- run_c8_gsea(setNames(gep_scores$X1, rownames(gep_scores)), species = "human", p_cutoff=1)
# res_gep2 <- run_c8_gsea(setNames(gep_scores$X2, rownames(gep_scores)), species = "human", p_cutoff=1)
# res_gep3 <- run_c8_gsea(setNames(gep_scores$X3, rownames(gep_scores)), species = "human", p_cutoff=1)
# res_gep4 <- run_c8_gsea(setNames(gep_scores$X4, rownames(gep_scores)), species = "human", p_cutoff=1)



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

p1 <- ggplot(gep_freqs, aes(
  x = Perturbation,
  y = freq,
  group = Composition_lvl2,
  color = Composition_lvl2
)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ GEP, scales = "free_y") +
  theme_bw(14) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = "Fraction of cells expressing GEP",
    color = "Group"
  ) + 
  scale_color_brewer(palette = "Set2")

p2 <- ggplot(gep_freqs, aes(
  x = Perturbation,
  y = freq,
  group = GEP,
  color = GEP
)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Composition_lvl2, scales = "free_y") +
  theme_bw(14) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = "Fraction of cells expressing GEP",
    color = "Group"
  ) + 
  scale_color_brewer(palette = "Set2")

p <- p1 + p2
ggsave(p, path="/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/figures", filename= "GEP_frequency.pdf", width=30, height=15, dpi = 700, units = "cm")

# U87: GEP4 goes up. According to GO GSEA of GEP 4, "negative regulation of gene expression, epigenetic" is enriched, which is consistent with mRNA-seq. However, GEP 2 (MES) does not go down.
# 