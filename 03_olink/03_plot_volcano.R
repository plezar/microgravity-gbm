library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)

da_res_list <- readRDS("03_olink/results/da_res_list.rds")

plot_volcano <- function(df) {
  df$gene <- rownames(df)
  
  # classification
  df <- df %>%
    mutate(
      negLogP = -log10(adj.P.Val),
      sig = case_when(
        adj.P.Val < 0.05 & logFC > 0 ~ "Up",
        adj.P.Val < 0.05 & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      sig = factor(sig, levels = c("Down", "NS", "Up"))
    )
  
  # pretty limits for plotting
  maxY <- max(df$negLogP, na.rm = TRUE)
  
  p <- ggplot(df, aes(x = logFC, y = negLogP)) +
    
    # significance cutoff line
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    
    # points
    geom_point(aes(color = sig), size = 2.2, alpha = 0.9) +
    
    # labels for significant genes
    geom_text_repel(
      data = subset(df, sig != "NS"),
      aes(label = gene),
      size = 4,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.size = 0.3,
      min.segment.length = 0
    ) +
    
    scale_color_manual(
      values = c(
        "Down" = "blue",
        "NS"   = "grey80",
        "Up"   = "red"
      )
    ) +
    
    labs(
      color = NULL,
      x = "logFC",
      y = expression(-log[10]("P"))
    ) +
    
    theme_bw(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      legend.position = "none",
      legend.justification = "left",
      legend.direction = "horizontal",
      legend.text = element_text(size = 16)
    )
  
  p
}

plot_list <- lapply(da_res_list, plot_volcano)
names(plot_list) <- names(da_res_list)


# give each plot a title
plot_list_titled <- Map(function(p, nm) p + ggtitle(nm),
                        plot_list,
                        names(plot_list))

# combine in grid (auto layout)
wrap_plots(plot_list_titled, ncol = 3)