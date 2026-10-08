#' Create a volcano plot for differential expression analysis
#'
#' This function creates a volcano plot to visualize differential expression results,
#' showing the relationship between fold change (x-axis) and statistical significance 
#' (y-axis). Points are colored based on significance thresholds, and the most 
#' significant genes can be labeled.
#'
#' @param results_df Data frame containing differential expression results
#' @param pval_colname Character string specifying the column name for nominal p-values (default: "pval")
#' @param padj_colname Character string specifying the column name for adjusted p-values (default: "padj")
#' @param logfc_colname Character string specifying the column name for log fold change (default: "logFC")
#' @param genesymbol_colname Character string specifying the column name for gene symbols (default: "GeneSymbol")
#' @param title Character string for plot title (default: NULL)
#' @param label_top_n Integer specifying number of top significant genes to label (default: 10)
#' @param plot_nominal_p Logical, if TRUE uses nominal p-values instead of adjusted p-values for y-axis (default: FALSE)
#' @param fc_threshold Numeric threshold for fold change significance (default: 1)
#' @param col_vec Named vector of colors for downregulated (-1), non-significant (0), and upregulated (1) genes
#' @param padj_threshold Numeric threshold for adjusted p-value significance (default: 0.05)
#'
#' @return A ggplot object representing the volcano plot
#' 
#' @details The function categorizes genes into three groups:
#' \itemize{
#'   \item Upregulated: log2FC > fc_threshold AND padj < padj_threshold (color: red)
#'   \item Downregulated: log2FC < -fc_threshold AND padj < padj_threshold (color: blue)
#'   \item Non-significant: all other genes (color: gray)
#' }
#' 
#' @examples
#' # Basic usage
#' volcano_plot(results_df)
#' 
#' # Custom column names and thresholds
#' volcano_plot(results_df, 
#'              logfc_colname = "coef",
#'              padj_colname = "FDR",
#'              fc_threshold = 0.5,
#'              padj_threshold = 0.01)
#'
#' @import ggplot2
#' @import dplyr
#' @import ggrepel
#' @import cowplot
#'
volcano_plot <- function(results_df,
                         pval_colname = "pval",
                         padj_colname = "padj",
                         logfc_colname = "logFC",
                         genesymbol_colname = "GeneSymbol",
                         title = NULL,
                         xlab = expression(log[2](FC)),
                         ylab = expression(-log[10](P)),
                         label_top_n = 10,
                         fc_threshold = 1,
                         ylim_min = 0,
                         genes_to_label = NULL,
                         ggrepel_max_overlaps = 5,
                         ggrepel_force = 50,
                         ylim_max,
                         xlim_min,
                         xlim_max,
                         col_vec = setNames(c("#377EB8", "#d3d3d3", "#E41A1C"), 
                                            c(-1, 0, 1)),
                         rect_col_vec = setNames(c("#377EB8", "#d3d3d3", "#E41A1C"), 
                                                 c(-1, 0, 1)),
                         padj_threshold = 0.05) {
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(cowplot)
  
  if (is.null(genes_to_label)) {
    genes_to_label <- results_df[[genesymbol_colname]]
  }
  
  # Prepare data with transformed variables and significance classification
  results_df_sub <- results_df %>%
    mutate(
      logPnominal = -log10(.data[[pval_colname]]),
      logP = -log10(.data[[padj_colname]]),
      log2FoldChange = .data[[logfc_colname]],
      gene_name = .data[[genesymbol_colname]],
      DE = case_when(
        log2FoldChange > fc_threshold & .data[[padj_colname]] < padj_threshold ~ 1,
        log2FoldChange < -fc_threshold & .data[[padj_colname]] < padj_threshold ~ -1,
        TRUE ~ 0
      )
    )
  
  # Use adjusted p-values (FDR)
  p <- ggplot(results_df_sub %>%
                arrange(desc(.data[[padj_colname]])),
              aes(x = log2FoldChange, y = logP, color = factor(DE))) +
    ## Right-upper (Up DEGs)
    annotate(
      "rect",
      xmin =  fc_threshold, xmax =  Inf,
      ymin = -log10(padj_threshold), ymax = Inf,
      fill  = rect_col_vec[[as.character(1)]],  # adjust if your "up" label isn't 1
      alpha = 0.50,
      colour = NA
    ) +
    ## Left-upper (Down DEGs)
    annotate(
      "rect",
      xmin = -Inf, xmax = -fc_threshold,
      ymin = -log10(padj_threshold), ymax = Inf,
      fill  = rect_col_vec[[as.character(-1)]], # adjust if your "down" label isn't -1
      alpha = 0.50,
      colour = NA
    ) +
    geom_point(size = 0.0001) +
    #geom_vline(xintercept = fc_threshold, lty = "dashed", linewidth = 0.5, alpha = 0.5) +
    #geom_vline(xintercept = -fc_threshold, lty = "dashed", linewidth = 0.5, alpha = 0.5) +
    #geom_hline(yintercept = -log10(padj_threshold), lty = "dashed", linewidth = 0.5, alpha = 0.5) +
    geom_text_repel(
      data = results_df_sub %>%
        filter(gene_name %in% genes_to_label) %>%
        filter(DE != 0) %>%
        arrange(desc(logP)) %>%
        head(label_top_n),
      aes(x = log2FoldChange, y = logP, label = gene_name),
      min.segment.length = 0,
      seed = 3,
      force = ggrepel_force,
      box.padding = 0.2,
      show.legend = FALSE,
      max.overlaps = ggrepel_max_overlaps,
      size = 5
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(ylim_min, ylim_max)
    ) +
    #coord_flip() +
    scale_color_manual(values = col_vec[names(col_vec) %in% unique(results_df_sub$DE)]) +
    theme_cowplot(14) +
    labs(x = xlab, y = ylab) +
    xlim(xlim_min, xlim_max) +
    ggtitle(title) +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.title = element_text(size = 12),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  
  return(p)
}