library(GSVA)
library(ggpubr)
library(rstatix)

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
  add_xy_position(x = "Group_gravity", fun = "mean_se")

p <- ggerrorplot(gsva_out_df, x = "Group_gravity", y = "score", 
            color = "Group_gravity",
            desc_stat = "mean_se",
            facet.by = c("Group_cell", "GEP"),
            add = "dotplot") +
  stat_pvalue_manual(stat.test, label = "p") +
  theme_bw(14)

ggsave(p, path="01_bulk-mrnaseq/results/figures", filename= "GEP_scores.pdf", width=30, height=15, dpi = 700, units = "cm")
