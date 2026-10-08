library(tidyverse)

u87_d <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87.csv", row.names = 1)
u87_d <- drop_na(u87_d)
write.table(data.frame(Acc = rownames(u87_d[u87_d$padj < 0.05, ])),
            "01_bulk-mrnaseq/results/homer_in/U87_DE_ALL.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(data.frame(Acc = rownames(u87_d[u87_d$padj < 0.05 & u87_d$log2FoldChange > 0, ])),
            "01_bulk-mrnaseq/results/homer_in/U87_DE_UP.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(data.frame(Acc = rownames(u87_d[u87_d$padj < 0.05 & u87_d$log2FoldChange < 0, ])),
            "01_bulk-mrnaseq/results/homer_in/U87_DE_DOWN.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

u87thp_d <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87THP1.csv", row.names = 1)
u87thp_d <- drop_na(u87thp_d)
write.table(data.frame(Acc = rownames(u87thp_d[u87thp_d$padj < 0.05, ])),
            "01_bulk-mrnaseq/results/homer_in/U87THP1_DE_ALL.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(data.frame(Acc = rownames(u87thp_d[u87thp_d$padj < 0.05 & u87thp_d$log2FoldChange > 0, ])),
            "01_bulk-mrnaseq/results/homer_in/U87THP1_DE_UP.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(data.frame(Acc = rownames(u87thp_d[u87thp_d$padj < 0.05 & u87thp_d$log2FoldChange < 0, ])),
            "01_bulk-mrnaseq/results/homer_in/U87THP1_DE_DOWN.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)