suppressPackageStartupMessages(library(metaseqR2))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cowplot))

count_mtx <- read.csv("01_bulk-mrnaseq/data/count_mtx.csv", row.names = 1)
meta_data <- read.csv("01_bulk-mrnaseq/data/meta_data.csv")
rownames(meta_data) <- meta_data$Sample
meta_data <- meta_data[colnames(count_mtx),]

## =====Remove duplicates======

rep1 <- "UuG1_b_S12_L001"
rep2 <- "UuG1_S1_L001"
new  <- "UuG1"   # name of collapsed sample

# ---- 1) sum counts ----
collapsed_counts <- count_mtx[, rep1] + count_mtx[, rep2]

# drop the two replicate columns and add the collapsed one
count_mtx2 <- count_mtx[, setdiff(colnames(count_mtx), c(rep1, rep2)), drop = FALSE]
count_mtx2 <- cbind(count_mtx2, collapsed_counts)
colnames(count_mtx2)[ncol(count_mtx2)] <- new


# ---- 2) collapse metadata ----
# Take one row as template (they should match for true tech reps)
md1 <- meta_data[rep1, , drop = FALSE]
md2 <- meta_data[rep2, , drop = FALSE]

meta_data2 <- meta_data[setdiff(rownames(meta_data), c(rep1, rep2)), , drop = FALSE]
md_new <- data.frame(Sample_Alias = "UuG1", Group_cell = "U87", Group_gravity = "uG", Sample = new)
rownames(md_new) <- new
meta_data2 <- rbind(meta_data2, md_new)

# ---- 3) final check: order meta_data to match count matrix ----
meta_data2 <- meta_data2[colnames(count_mtx2), , drop = FALSE]

# outputs:
# count_mtx2, meta_data2

## =====DESeq2=====

set.seed(42)
#count_mtx_downsampled <- downsampleCounts(count_mtx2)
dds <- DESeqDataSetFromMatrix(countData = count_mtx2,
                              colData = meta_data2,
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


### ===== Interaction =====

res_int <- results(dds, name = "Group_gravityuG.Group_cellU87THP")
res_int <- as.data.frame(res_int) %>% arrange(padj)
write.csv(res_int, "01_bulk-mrnaseq/results/deseq2_interaction.csv")




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

deg_counts_long$comparison <- factor(deg_counts_long$comparison, levels = c("U87", "U87THP", "uG", "KSC"))
levels(deg_counts_long$comparison) <- c("U87 (uG vs KSC)", "U87THP (uG vs KSC)", "uG (U87THP vs U87)", "KSC (U87THP vs U87)")
deg_counts_long$y <- deg_counts_long$count
deg_counts_long$y <- if_else(deg_counts_long$y < 9, deg_counts_long$y + 25, deg_counts_long$y + 1)

p <- ggplot(deg_counts_long, aes(x = comparison, y = count, fill = direction)) +
  geom_bar(stat = "identity", color = "black") +
  
  geom_text(
    aes(
      label = count,
      color = direction,
      y = y
    ),
    position = position_stack(vjust = 0.5),
    angle = 90,
    size = 3
  ) +
  
  scale_fill_manual(values = c("up" = "#f6d8cb", "down" = "#d6e6f2")) +
  scale_color_manual(values = c("up" = "#7b0f0f", "down" = "#0b3c6f")) +
  
  coord_flip() +
  scale_x_discrete(position = "top") +
  
  theme_cowplot(14) +
  labs(
    x = NULL,
    y = "DEG count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    )
  )

ggsave(p, path="01_bulk-mrnaseq/results/figures", filename= "DEG_counts.pdf", width=12, height=5, dpi = 700, units = "cm")


### ==== upsetr plot ======

# get gene sets for each comparison
DEG_sets <- list(
  "U87 (uG vs KSC)"    = rownames(res_U87)[which(res_U87$padj < 0.05)],
  "U87THP (uG vs KSC)" = rownames(res_U87THP)[which(res_U87THP$padj < 0.05)],
  "uG (U87THP vs U87)"     = rownames(res_uG)[which(res_uG$padj < 0.05)],
  "KSC (U87THP vs U87)"    = rownames(res_KSC)[which(res_KSC$padj < 0.05)]
)
pdf("01_bulk-mrnaseq/results/figures/deg_upset_plot.pdf", width = 3.5, height = 4)
upset(
  fromList(DEG_sets),
  nsets = 4,
  #mb.ratio = c(0.8, 0.5),
  order.by = "freq",
  set_size.show = F,
  sets = names(DEG_sets)
)
dev.off()
