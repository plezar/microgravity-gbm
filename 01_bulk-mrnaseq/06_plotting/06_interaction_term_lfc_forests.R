library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(cowplot)

genes <- c("FAM50A", "PRSS35", "EPHA7", "LINC01705")

df_u87 <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87.csv", row.names = 1)
df_u87$gene <- rownames(df_u87)
df_u87$contrast <- "U87"

df_u87thp1 <- read.csv("01_bulk-mrnaseq/results/deseq2_uG_vs_KSC_U87THP1.csv", row.names = 1)
df_u87thp1$gene <- rownames(df_u87thp1)
df_u87thp1$contrast <- "U87THP1"

df <- bind_rows(df_u87, df_u87thp1) %>%
  filter(gene %in% genes) %>%
  mutate(
    gene = factor(gene, levels = genes),
    contrast = factor(contrast, levels = c("U87", "U87THP1"))
  )

levels(df$contrast) <- c("GBM", "+Mono")

# LFC +/- SE plot (interaction = flips/differences across contrasts)
p <- ggplot(df, aes(x = contrast, y = log2FoldChange, group = gene)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(
    aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE),
    width = 0.15, linewidth = 0.7
  ) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  facet_wrap(~ gene, nrow = 1) +
  labs(
    x = NULL,
    y = bquote(log[2](FC) %+-% SE)
  ) +
  theme_cowplot() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggsave(p, path="01_bulk-mrnaseq/results/figures", filename= "interaction_dotplot_sem.pdf", width=10, height=6.5, dpi = 700, units = "cm")
