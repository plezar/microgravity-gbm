suppressPackageStartupMessages(library(metaseqR2))
suppressPackageStartupMessages(library(DESeq2))

count_mtx <- read.csv("01_bulk-mrnaseq/data/count_mtx.csv", row.names = 1)
meta_data <- read.csv("01_bulk-mrnaseq/data/meta_data.csv")
rownames(meta_data) <- meta_data$Sample
meta_data <- meta_data[colnames(count_mtx),]

## =====DESeq2=====

set.seed(42)
count_mtx_downsampled <- downsampleCounts(count_mtx)
dds <- DESeqDataSetFromMatrix(countData = count_mtx_downsampled,
                              colData = meta_data,
                              design =  ~ Group_gravity + Group_cell + Group_gravity:Group_cell)

# pre-filtering
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds)
saveRDS(dds, "01_bulk-mrnaseq/results/dds.rds")

# resultsNames(dds)
# Intercept: expression in the baseline (U87, KSC)
# Group_gravity_uG_vs_KSC: Microgravity effects in ref. (U87)
# Group_cell_U87THP_vs_U87: THP1 effects in ref. (KSC)
# Group_gravityuG.Group_cellU87THP: change in LFC when going from ref (U87) to U87+THP1


### =====Microgravity effect=====

res_U87 <- results(dds, name = "Group_gravity_uG_vs_KSC")
res_U87 <- as.data.frame(res_U87) %>% arrange(padj)
write.csv(res_U87, "01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87.csv")


res_U87THP <- results(dds, contrast = list(
  c("Group_gravity_uG_vs_KSC", "Group_gravityuG.Group_cellU87THP")
))
res_U87THP <- as.data.frame(res_U87THP) %>% arrange(padj)
write.csv(res_U87THP, "01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87THP1.csv")


### =====Co-culture effect=====

res_KSC <- results(dds, name = "Group_cell_U87THP_vs_U87")
res_KSC <- as.data.frame(res_KSC) %>% arrange(padj)
write.csv(res_KSC, "01_bulk-mrnaseq/results/deseq2_U87THP_vs_U87_KSC.csv")


res_uG <- results(dds, contrast = list(
  c("Group_cell_U87THP_vs_U87", "Group_gravityuG.Group_cellU87THP")
))
res_uG <- as.data.frame(res_uG) %>% arrange(padj)
write.csv(res_uG, "01_bulk-mrnaseq/results/deseq2_U87THP_vs_U87_uG.csv")



## ====count DEGs=====

library(tidyverse)
library(UpSetR)

# helper function
count_deg <- function(res) {
  res %>% 
    filter(!is.na(padj), padj < 0.05) %>%
    summarise(
      up = sum(log2FoldChange > 0, na.rm = TRUE),
      down = sum(log2FoldChange < 0, na.rm = TRUE)
    )
}

deg_counts <- bind_rows(
  U87      = count_deg(res_U87),
  U87THP   = count_deg(res_U87THP),
  KSC      = count_deg(res_KSC),
  uG       = count_deg(res_uG),
  .id = "comparison"
)


### ==== bar plot ======

deg_counts_long <- deg_counts %>%
  pivot_longer(cols = c(up, down),
               names_to = "direction",
               values_to = "count")

p <- ggplot(deg_counts_long, aes(x = comparison, y = count, fill = direction)) +
  geom_bar(stat = "identity") +   # <-- stacked by default
  scale_fill_manual(values = c("up" = "firebrick", "down" = "steelblue")) +
  theme_bw(14) +
  labs(x = "Comparison",
       y = "DEG count")
ggsave(p, path="01_bulk-mrnaseq/results/figures", filename= "DEG_counts.pdf", width=10, height=5, dpi = 700, units = "cm")


### ==== upsetr plot ======

# get gene sets for each comparison
DEG_sets <- list(
  U87    = rownames(res_U87)[which(res_U87$padj < 0.05)],
  U87THP = rownames(res_U87THP)[which(res_U87THP$padj < 0.05)],
  KSC    = rownames(res_KSC)[which(res_KSC$padj < 0.05)],
  uG     = rownames(res_uG)[which(res_uG$padj < 0.05)]
)

pdf("01_bulk-mrnaseq/results/figures/deg_upset_plot.pdf", width = 4, height = 3)
upset(
  fromList(DEG_sets),
  nsets = 4,
  order.by = "freq",
  sets = names(DEG_sets)
)
dev.off()
