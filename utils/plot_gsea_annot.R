wrap_at_underscore <- function(x, width = 25, sep = "_") {
  x <- as.character(x)
  
  vapply(x, function(s) {
    parts <- strsplit(s, sep, fixed = TRUE)[[1]]
    if (length(parts) <= 1) return(s)
    
    lines <- character(0)
    cur <- parts[1]
    
    for (i in 2:length(parts)) {
      candidate <- paste0(cur, sep, parts[i])
      
      if (nchar(candidate) <= width) {
        cur <- candidate
      } else {
        lines <- c(lines, cur)
        cur <- parts[i]
      }
    }
    
    lines <- c(lines, cur)
    paste(lines, collapse = "\n")
  }, character(1))
}

plot_gsea_annot <- function(gsea_cp,
                            path,
                            color = "#7b0f0f",
                            subplots = 1:2,
                            title = path,
                            label_size = 5,
                            label_fill = "white",
                            label_alpha = 0.75,
                            label_color = "black",
                            hjust = 1.05,
                            vjust = 1.15,
                            wrap_width = 25,
                            show_legend = FALSE,
                            show_x_axis_elements = TRUE) {
  stopifnot(!missing(gsea_cp), !missing(path))
  
  # pull results table
  gsea_cp_df <- as.data.frame(gsea_cp@result)
  
  if (!path %in% gsea_cp_df$ID) {
    stop("Pathway not found in gsea_cp results: ", path)
  }
  
  lab_df <- gsea_cp_df %>%
    filter(ID == path) %>%
    mutate(path = wrap_at_underscore(path, width = wrap_width)) %>%
    transmute(
      x = Inf, y = Inf,
      label = sprintf("%s\nNES = %.2f\nP = %.3g", path, NES, p.adjust)
    )
  
  p <- enrichplot::gseaplot2(
    gsea_cp,
    geneSetID = path,
    title = title,
    color = color,
    subplots = subplots
  )
  
  # annotate running ES panel (first subplot)
  p[[1]] <- p[[1]] +
    theme(
      axis.line = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    ylab("ES") +
    geom_label(
      data = lab_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = hjust, vjust = vjust,
      size = label_size,
      label.size = 0,
      fill = label_fill,
      alpha = label_alpha,
      color = label_color
    ) +
    ggtitle(NULL)
  
  # optional cleanup of the hits panel (second subplot), if present
  if (length(p) >= 2) {
    p[[2]] <- p[[2]] + theme(axis.line.x = element_blank())
  }
  
  if (!show_x_axis_elements) {
    p[[2]] <- p[[2]] + theme(axis.text.x = element_blank(),
                             axis.ticks.x = element_blank())
  }
  
  p
}