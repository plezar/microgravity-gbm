library(GSVA)
library(ggpubr)
library(rstatix)
library(DESeq2)
library(tidyverse)
library(cowplot)
source("utils/utils.R")

top_genes <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/xenium_cNMF_k4_topgenes.csv", row.names = 1)
top_genes <- head(top_genes, 50)
top_genes <- as.list(top_genes)
names(top_genes) <- paste0("GEP ", 1:4)

dds <- readRDS("01_bulk-mrnaseq/results/dds.rds")

norm_mat <- assay(vst(dds))
param    <- ssgseaParam(norm_mat, top_genes)
gsva_out <- gsva(param)
gsva_out <- t(scale(t(gsva_out)))

gsva_out_df <- t(gsva_out) %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(!Sample, names_to = "GEP", values_to = "score") %>%
  left_join(as.data.frame(colData(dds)))


stat.test <- gsva_out_df %>%
  group_by(GEP, Group_cell) %>%
  rstatix::wilcox_test(score ~ Group_gravity) %>%
  add_significance("p") %>%
  add_xy_position(x = "Group_gravity", fun = "mean_se") %>%
  adjust_pvalue()

levels(gsva_out_df$Group_cell) <- c("GBM", "+Mono")

p <- ggerrorplot(gsva_out_df, x = "Group_gravity", y = "score", 
            color = "Group_gravity",
            desc_stat = "mean_se",
            #facet.by = c("Group_cell", "GEP"),
            add = "dotplot") +
  stat_pvalue_manual(stat.test, label = "p.adj", hide.ns = T) +
  facet_grid(Group_cell ~ GEP) +
  theme_cowplot(14) +
  scale_color_manual(values = COLPAL_UG_KSC) +
  labs(x = NULL, y = "ssGSEA score") +
  theme(
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(p, path="02_st-xenium/figures", filename= "GEP_scores.pdf", width=11*1.2, height=8*1.2, dpi = 700, units = "cm")
