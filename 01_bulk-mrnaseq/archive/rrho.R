library(RRHO2)
library(RColorBrewer)

run_RRHO2 <- function(df1, df2, label1 = "Sample1", label2 = "Sample2", rank_col) {
  
  df1$symbol <- rownames(df1)
  df2$symbol <- rownames(df2)
  
  df1 <- df1[, c("symbol", rank_col)]
  df2 <- df2[, c("symbol", rank_col)]
  
  shared_genes <- intersect(df1$symbol, df2$symbol)
  
  df1_sub <- df1[df1$symbol %in% shared_genes, ]
  df2_sub <- df2[df2$symbol %in% shared_genes, ]
  
  df1_sub <- df1_sub[match(shared_genes, df1_sub$symbol), ]
  df2_sub <- df2_sub[match(shared_genes, df2_sub$symbol), ]
  
  rr_up <- RRHO2_initialize(df1_sub,
                            df2_sub,
                            labels = c(label1, label2),
                            stepsize = nrow(df1_sub) / 100,
                            method = "hyper")
  
  return(rr_up)
}

df_u87 <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87.csv", row.names = 1)
df_u87thp1 <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87THP1.csv", row.names = 1)

spectral_continuous <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(256)

rr <- run_RRHO2(df_u87, df_u87thp1, "U87", "U87+Mono", "log2FoldChange")

pdf("Figure_3/figures/S3B_iN.pdf", width = 6, height = 5.5)
RRHO2_heatmap(rr, colorGradient = spectral_continuous)
dev.off()


# correlation

df_u87$symbol <- rownames(df_u87)
df_u87thp1$symbol <- rownames(df_u87thp1)

df_u87 <- df_u87[, c("symbol", "log2FoldChange")]
df_u87thp1 <- df_u87thp1[, c("symbol", "log2FoldChange")]

shared_genes <- intersect(df_u87$symbol, df_u87thp1$symbol)

df1_sub <- df_u87[df_u87$symbol %in% shared_genes, ]
df2_sub <- df_u87thp1[df_u87thp1$symbol %in% shared_genes, ]

df1_sub <- df1_sub[match(shared_genes, df1_sub$symbol), ]
df2_sub <- df2_sub[match(shared_genes, df2_sub$symbol), ]

cor(df1_sub$log2FoldChange, df2_sub$log2FoldChange, method = "spearman")
